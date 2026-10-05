"""Turn raw OCHEM property exports into the tidy file read by Talktorial T033.

OCHEM (https://ochem.eu) requires an account, so its data cannot be downloaded
from inside a teaching notebook. Instead, the exports were downloaded once,
preprocessed with this script, and the (small) result was committed to the
repository -- the same pattern Talktorial T032 uses for its ClustalO alignment.

Reading an OCHEM export
-----------------------
An export holds two value columns per property:

    "<Property>"                        + "UNIT {<Property>}"
    "<Property> {measured, converted}"  + "UNIT {<Property>}.1"

The first is the value exactly as the source reported it, so its unit column is a
mixture ('?C', 'K', 'mm Hg', 'log(mmHg)', 'mg/L', ...). The second is the value
after OCHEM's on-the-fly conversion, and it is homogeneous. **We use only the
converted column and its own unit column**, and assert that the unit is the one
we expect before using the numbers.

What this script does
---------------------
1. Reads every raw export in ``data/raw`` and groups them by property; a property
   split over several files (OCHEM caps the size of a single download) is
   concatenated and de-duplicated.
2. Checks that the converted unit is the expected one, and applies the derived
   transform the notebook needs (see ``PROPERTIES`` below).
3. Standardises every structure with RDKit: largest organic fragment, neutralised,
   isotopes and stereochemistry removed; the InChIKey is the compound identity.
4. Drops records that are not usable (unparsable SMILES, missing values,
   inorganics, very large molecules, values outside a physically sensible range).
5. Aggregates the remaining records per (compound, property) into a median value,
   keeping the number of records and their spread so that the notebook can talk
   about experimental uncertainty.
6. Writes a long-format, gzip-compressed CSV to ``data/ochem_properties.csv.gz``.

Usage
-----
    python scripts/prepare_ochem_data.py                 # use data/raw
    python scripts/prepare_ochem_data.py --raw-dir <dir> # somewhere else
    python scripts/prepare_ochem_data.py --inspect       # only report what was found
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

# --------------------------------------------------------------------------- #
# What we expect to find in data/raw.                                          #
#                                                                              #
# ``file``        substring matched against the raw file names, after lowercasing
#                 and stripping spaces/underscores -- so one entry picks up both
#                 ``meltingpoint.csv`` and ``meltingpoint_2.csv``.             #
# ``source_unit`` the unit OCHEM's converted column is in. Checked, not assumed.#
# ``transform``   what we do to get from ``source_unit`` to ``unit``.           #
# ``limits``      sanity range in the final unit; values outside it are almost  #
#                 always unit mix-ups or typos in the underlying literature.    #
# --------------------------------------------------------------------------- #
PROPERTIES = {
    "melting_point": {
        "file": "meltingpoint",
        "ochem_name": "Melting Point",
        "source_unit": "Celsius",
        "transform": None,
        "unit": "°C",
        "label": "Melting point",
        "limits": (-200.0, 600.0),
    },
    "logp": {
        "file": "logpow",
        "ochem_name": "logPow",
        "source_unit": "Log unit",
        "transform": None,
        "unit": "log unit",
        "label": "logP (octanol/water)",
        "limits": (-6.0, 12.0),
    },
    "logs": {
        "file": "watersolubility",
        "ochem_name": "Water solubility",
        # OCHEM converts water solubility to *linear* mol/L; the modelling
        # literature (and the General Solubility Equation) works in log units, so
        # we take the base-10 logarithm here. Round-tripping the records whose
        # source unit was already log(mol/L) agrees to ~3e-5 log units, so this
        # costs no precision.
        "source_unit": "mol/L",
        "transform": "log10",
        "unit": "log(mol/L)",
        "label": "logS (water solubility)",
        "limits": (-12.0, 2.0),
    },
    "boiling_point": {
        "file": "boilingpoint",
        "ochem_name": "Boiling Point",
        "source_unit": "Celsius",
        "transform": None,
        "unit": "°C",
        "label": "Boiling point",
        "limits": (-150.0, 700.0),
    },
    "vapour_pressure": {
        "file": "vaporpressure",
        "ochem_name": "Vapor Pressure",
        # OCHEM converts vapour pressure to log(Pa). We keep that SI unit rather
        # than shifting to log(mmHg), so no conversion is applied here.
        "source_unit": "log(Pa)",
        "transform": None,
        "unit": "log(Pa)",
        "label": "log vapour pressure",
        "limits": (-14.0, 8.0),
    },
}

MAX_HEAVY_ATOMS = 70
MIN_HEAVY_ATOMS = 3

# Column names OCHEM has used for the identifiers we keep, lowercased.
SMILES_ALIASES = ("smiles", "canonical smiles", "structure")
CASRN_ALIASES = ("casrn", "cas", "cas number")
ARTICLE_ALIASES = ("articleid", "article id", "article")


# --------------------------------------------------------------------------- #
# Column discovery                                                             #
# --------------------------------------------------------------------------- #
def _find_column(columns, aliases):
    lowered = {str(c).strip().lower(): c for c in columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    for alias in sorted(aliases, key=len, reverse=True):
        for low, original in lowered.items():
            if alias in low:
                return original
    return None


def _converted_columns(columns, ochem_name):
    """Locate the '{measured, converted}' value column and *its* unit column.

    pandas disambiguates the two identically named ``UNIT {<Property>}`` headers
    by suffixing the second one with ``.1``; the converted value belongs to the
    last of them.
    """
    value_column = next(
        (c for c in columns if str(c).strip() == f"{ochem_name} {{measured, converted}}"),
        None,
    )
    unit_pattern = re.compile(rf"UNIT \{{{re.escape(ochem_name)}\}}(\.\d+)?$")
    unit_columns = [c for c in columns if unit_pattern.fullmatch(str(c).strip())]
    unit_column = unit_columns[-1] if unit_columns else None
    return value_column, unit_column


def read_raw(path):
    """Read one OCHEM export, whatever tabular format it was saved in."""
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path)
    for sep in (",", "\t", ";"):
        frame = pd.read_csv(path, sep=sep, low_memory=False)
        if frame.shape[1] > 1:
            return frame
    raise ValueError(f"could not work out the delimiter of {path}")


# --------------------------------------------------------------------------- #
# Structure standardisation                                                    #
# --------------------------------------------------------------------------- #
class Standardiser:
    """Cache-backed SMILES -> (canonical SMILES, InChIKey) standardisation."""

    def __init__(self):
        self._chooser = rdMolStandardize.LargestFragmentChooser()
        self._uncharger = rdMolStandardize.Uncharger()
        self._cache: dict[str, tuple[str, str] | None] = {}

    def __call__(self, smiles):
        if smiles in self._cache:
            return self._cache[smiles]
        result = self._standardise(smiles)
        self._cache[smiles] = result
        return result

    def _standardise(self, smiles):
        if not isinstance(smiles, str) or not smiles.strip():
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        try:
            mol = rdMolStandardize.Cleanup(mol)
            mol = self._chooser.choose(mol)
            mol = self._uncharger.uncharge(mol)
        except Exception:
            return None
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        if not MIN_HEAVY_ATOMS <= mol.GetNumHeavyAtoms() <= MAX_HEAVY_ATOMS:
            return None
        if not any(atom.GetSymbol() == "C" for atom in mol.GetAtoms()):
            return None  # keep the set organic
        # no isotopes, no stereochemistry: OCHEM measurements are rarely resolved
        # at that level and it would needlessly fragment the compound identity
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        Chem.RemoveStereochemistry(mol)
        try:
            inchikey = Chem.MolToInchiKey(mol)
        except Exception:
            return None
        if not inchikey:
            return None
        smiles_out = Chem.MolToSmiles(mol)
        # Insist that the result parses again. Stripping the metal off an
        # organometallic (ferrocenes, mostly) can leave an aromatic
        # cyclopentadienyl ring that RDKit writes but refuses to read back, and
        # a SMILES in the output file that nothing can parse is a trap for every
        # downstream user.
        if Chem.MolFromSmiles(smiles_out) is None:
            return None
        return smiles_out, inchikey


# --------------------------------------------------------------------------- #
# Per-property pipeline                                                        #
# --------------------------------------------------------------------------- #
def locate_files(raw_dir):
    """Map each property onto every raw file that seems to hold it."""
    files = [p for p in sorted(raw_dir.iterdir())
             if p.suffix.lower() in (".csv", ".tsv", ".txt", ".xlsx", ".xls")]
    found = {}
    for key, spec in PROPERTIES.items():
        def normalise(name):
            return name.lower().replace(" ", "").replace("_", "").replace("-", "")
        matches = [p for p in files if spec["file"] in normalise(p.name)]
        if matches:
            found[key] = matches
    return found, files


def collect_records(key, paths, verbose=True):
    """Read every file of one property and return a tidy (smiles, value) table."""
    spec = PROPERTIES[key]
    pieces = []

    for path in paths:
        frame = read_raw(path)
        value_column, unit_column = _converted_columns(frame.columns, spec["ochem_name"])
        smiles_column = _find_column(frame.columns, SMILES_ALIASES)
        casrn_column = _find_column(frame.columns, CASRN_ALIASES)
        article_column = _find_column(frame.columns, ARTICLE_ALIASES)

        if verbose:
            print(f"\n[{key}] {path.name}: {len(frame):,} rows, {frame.shape[1]} columns")
            print(f"    value -> {value_column!r}")
            print(f"    unit  -> {unit_column!r}")
            print(f"    smiles-> {smiles_column!r}")

        if value_column is None or smiles_column is None:
            raise SystemExit(
                f"[{key}] {path.name}: could not find the "
                f"'{spec['ochem_name']} {{measured, converted}}' and/or SMILES column.\n"
                f"Columns present: {list(frame.columns)}"
            )

        # The unit of the converted column is checked, never assumed: OCHEM
        # converts on the fly and a different export setting would silently
        # change the numbers.
        units = set(frame[unit_column].dropna().astype(str).str.strip()) if unit_column else set()
        if verbose:
            print(f"    converted unit(s): {sorted(units) or 'none reported'}")
        if units and units != {spec["source_unit"]}:
            raise SystemExit(
                f"[{key}] {path.name}: expected the converted column to be in "
                f"'{spec['source_unit']}', found {sorted(units)}. Either re-export "
                f"this property in {spec['source_unit']}, or update "
                f"PROPERTIES['{key}']['source_unit'] and ['transform']."
            )

        pieces.append(pd.DataFrame({
            "smiles_raw": frame[smiles_column],
            "value": pd.to_numeric(frame[value_column], errors="coerce"),
            "casrn": frame[casrn_column] if casrn_column else pd.NA,
            "article": frame[article_column] if article_column else pd.NA,
        }))

    records = pd.concat(pieces, ignore_index=True)
    if len(paths) > 1:
        before = len(records)
        records = records.drop_duplicates(subset=["smiles_raw", "value", "article"])
        if verbose:
            print(f"\n[{key}] {before - len(records):,} duplicate records removed "
                  f"where the split downloads overlapped")
    return records


def process_property(key, paths, standardiser, verbose=True):
    spec = PROPERTIES[key]
    records = collect_records(key, paths, verbose=verbose)
    n_start = len(records)

    records = records[records["value"].notna()]

    if spec["transform"] == "log10":
        positive = records["value"] > 0
        if verbose:
            print(f"[{key}] {(~positive).sum():,} non-positive values dropped before log10")
        records = records[positive].copy()
        records["value"] = np.log10(records["value"])
    elif spec["transform"] is not None:
        raise ValueError(f"unknown transform {spec['transform']!r}")

    low, high = spec["limits"]
    in_range = records["value"].between(low, high)
    if verbose:
        print(f"[{key}] {(~in_range).sum():,} values outside {low} to {high} "
              f"{spec['unit']} dropped")
    records = records[in_range]

    standardised = records["smiles_raw"].map(standardiser)
    keep = standardised.notna()
    records = records[keep].copy()
    records["smiles"] = [s[0] for s in standardised[keep]]
    records["inchikey"] = [s[1] for s in standardised[keep]]

    aggregated = (
        records.groupby("inchikey")
        .agg(
            smiles=("smiles", "first"),
            casrn=("casrn", "first"),
            value=("value", "median"),
            n_records=("value", "size"),
            value_std=("value", "std"),
        )
        .reset_index()
    )
    aggregated["property"] = key
    aggregated["unit"] = spec["unit"]

    if verbose:
        print(f"[{key}] kept {len(records):,} of {n_start:,} records "
              f"-> {len(aggregated):,} unique compounds")
    return aggregated


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=DATA / "raw")
    parser.add_argument("--out", type=Path, default=DATA / "ochem_properties.csv.gz")
    parser.add_argument("--min-tasks", type=int, default=1,
                        help="only keep compounds measured for at least this many properties")
    parser.add_argument("--max-single-task", type=int, default=None,
                        help="cap how many compounds carrying only one measurement are "
                             "kept, to keep the committed file small")
    parser.add_argument("--seed", type=int, default=22)
    parser.add_argument("--inspect", action="store_true",
                        help="report the files, columns and units found, then stop")
    args = parser.parse_args(argv)

    if not args.raw_dir.is_dir():
        raise SystemExit(f"raw directory {args.raw_dir} does not exist")

    found, all_files = locate_files(args.raw_dir)
    missing = sorted(set(PROPERTIES) - set(found))
    print(f"raw directory: {args.raw_dir}")
    print(f"files present: {[p.name for p in all_files]}")
    for key, paths in found.items():
        print(f"  {key:16s} <- {[p.name for p in paths]}")
    if missing:
        print(f"not found    : {missing} (these properties will be skipped)")
    if not found:
        raise SystemExit("no raw property exports found")

    if args.inspect:
        standardiser = None
        for key, paths in found.items():
            collect_records(key, paths, verbose=True)
        return

    standardiser = Standardiser()
    tables = [process_property(key, paths, standardiser) for key, paths in found.items()]
    tidy = pd.concat(tables, ignore_index=True)

    counts = tidy.groupby("inchikey")["property"].nunique()
    if args.min_tasks > 1:
        tidy = tidy[tidy["inchikey"].isin(counts[counts >= args.min_tasks].index)]

    if args.max_single_task is not None:
        single = counts[counts == 1].index
        if len(single) > args.max_single_task:
            generator = np.random.default_rng(args.seed)
            drop = generator.choice(single.to_numpy(),
                                    size=len(single) - args.max_single_task,
                                    replace=False)
            tidy = tidy[~tidy["inchikey"].isin(drop)]
            print(f"\nkept a random {args.max_single_task:,} of the {len(single):,} "
                  f"compounds that carry only one measurement")

    tidy = tidy[["inchikey", "smiles", "casrn", "property", "value", "unit",
                 "n_records", "value_std"]]
    tidy = tidy.sort_values(["property", "inchikey"]).reset_index(drop=True)
    tidy["value"] = tidy["value"].round(3)
    tidy["value_std"] = tidy["value_std"].round(3)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tidy.to_csv(args.out, index=False, compression="gzip")

    print("\n" + "=" * 62)
    print(f"wrote {args.out}  ({args.out.stat().st_size / 1e6:.1f} MB)")
    print(f"{len(tidy):,} compound-property pairs, "
          f"{tidy['inchikey'].nunique():,} unique compounds")
    print("\nrecords per property:")
    print(tidy.groupby(["property", "unit"]).size().to_string())
    print("\ncompounds by number of properties measured:")
    print(tidy.groupby("inchikey")["property"].nunique()
          .value_counts().sort_index().to_string())


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Set, Tuple


SUPPORTED_EXTENSIONS: frozenset = frozenset({'.xlsx', '.xls', '.csv', '.tsv', '.tab'})


def load_file(filepath: str) -> pd.DataFrame:
    path = Path(filepath)
    ext  = path.suffix.lower()

    if ext not in SUPPORTED_EXTENSIONS:
        sys.exit(
            f"ERROR: Unsupported file extension '{ext}' — {filepath}\n"
            f"Supported formats: {', '.join(sorted(SUPPORTED_EXTENSIONS))}"
        )

    if not path.is_file():
        sys.exit(f"ERROR: File not found: {filepath}")

    try:
        if ext in ('.xlsx', '.xls'):
            df = pd.read_excel(filepath, header=0, dtype=str)
        elif ext == '.csv':
            df = pd.read_csv(filepath, header=0, dtype=str)
        else:  # .tsv / .tab
            df = pd.read_csv(filepath, sep='\t', header=0, dtype=str)
    except Exception as exc:
        sys.exit(f"ERROR: Could not read '{filepath}': {exc}")

    if df.empty:
        sys.exit(f"ERROR: File is empty or has no data rows: {filepath}")

    return df


def merge_files(
    files: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    warnings: List[str] = []

    if not files:
        sys.exit("ERROR: No input files provided.")

    # Load and normalise the Sample_ID column
    loaded: List[Tuple[str, pd.DataFrame]] = []
    seen_stems: Dict[str, str] = {}

    for filepath in files:
        df   = load_file(filepath)
        stem = Path(filepath).stem

        if stem in seen_stems:
            warnings.append(
                f"WARNING: Duplicate filename stem '{stem}' "
                f"('{filepath}' vs '{seen_stems[stem]}'). "
                "Duplicate-column suffixes may collide."
            )
        seen_stems[stem] = filepath

        # Rename first column → 'Sample_ID'
        original_col = str(df.columns[0])
        if original_col != 'Sample_ID':

            if 'Sample_ID' in df.columns:
                warnings.append(
                    f"WARNING: '{filepath}' has a non-first column named "
                    f"'Sample_ID'. It will be renamed 'Sample_ID_{stem}' "
                    "to avoid conflict with the join key."
                )
                df = df.rename(columns={'Sample_ID': f'Sample_ID_{stem}'})
            df = df.rename(columns={original_col: 'Sample_ID'})

        # Strip accidental whitespace from Sample_ID values.
        df['Sample_ID'] = df['Sample_ID'].str.strip()

        loaded.append((stem, df))

    # Find column names shared by more than one file 
    col_count: Dict[str, int] = {}
    for _, df in loaded:
        for col in df.columns.tolist():
            if col != 'Sample_ID':
                col_count[col] = col_count.get(col, 0) + 1

    dup_cols: Set[str] = {col for col, n in col_count.items() if n > 1}

    if dup_cols:
        warnings.append(
            "NOTE: Duplicate column header(s) across files — "
            "appending filename stem as suffix: "
            + ', '.join(sorted(dup_cols))
        )

    # Rename duplicate columns symmetrically in every affected file 
    resolved: List[pd.DataFrame] = []
    for stem, df in loaded:
        rename_map = {
            col: f"{col}_{stem}"
            for col in df.columns.tolist()
            if col != 'Sample_ID' and col in dup_cols
        }
        resolved.append(df.rename(columns=rename_map) if rename_map else df)

    # Outer-merge all frames progressively on Sample_ID
    merged = resolved[0]
    for df in resolved[1:]:
        merged = pd.merge(merged, df, on='Sample_ID', how='outer')

    # Replace every NaN with the sentinel string 'NA' 
    merged = merged.fillna('NA')

    return merged, warnings


def main() -> None:
    parser = argparse.ArgumentParser(
        prog='meta_merge',
        description=(
            "Merge metadata and genomic result files into a single CSV table "
            "using the first column of each file as the Sample_ID join key."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Supported formats: .xlsx  .xls  .csv  .tsv  .tab\n\n"
            "Example:\n"
            "  meta_merge.py metadata.xlsx quast_summary.tsv checkm2.csv \\\n"
            "                -o meta_merged.csv"
        ),
    )
    parser.add_argument(
        'files',
        nargs='+',
        help="Input files to merge (accepts mixed formats).",
    )
    parser.add_argument(
        '-o', '--output',
        default='meta_merged.csv',
        metavar='FILE',
        help='Output filename (default: meta_merged.csv)',
    )
    args = parser.parse_args()

    print(f"[meta_merge] Merging {len(args.files)} file(s):")
    for f in args.files:
        print(f"    {f}")

    merged, warnings = merge_files(args.files)

    for msg in warnings:
        print(msg, file=sys.stderr)

    merged.to_csv(args.output, index=False)

    print(
        f"[meta_merge] Done — {len(merged)} sample(s), "
        f"{len(merged.columns)} column(s) written to: {args.output}"
    )


if __name__ == '__main__':
    main()
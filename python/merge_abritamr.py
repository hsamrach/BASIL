#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from collections import Counter
from pathlib import Path

import pandas as pd

DEFAULT_FILENAME = "amrfinder.out"
DEFAULT_OUTPUT = "merged_abritamr.tsv"
DEFAULT_MINCOV_GENES = 0.0
DEFAULT_MINID_GENES = 0.0

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="merge_abritamr.py",
        description=(
            "script for merging abritamr output and filtering its coverages and identities"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        type=Path,
        metavar="DIR",
        help="Input directory containing per-sample subdirectories with abritamr output files (amrfinder.out).",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        metavar="FILE",
        help="Path of the merged output TSV.",
    )
    parser.add_argument(
        "--filename",
        default=DEFAULT_FILENAME,
        metavar="NAME",
        help="Filename to search for within each sample directory.",
    )
    parser.add_argument(
        "--mincov_genes",
        type=float,
        default=DEFAULT_MINCOV_GENES,
        metavar="FLOAT",
        help="Minimum filtering coverages",
    )
    parser.add_argument(
        "--minid_genes",
        type=float,
        default=DEFAULT_MINID_GENES,
        metavar="FLOAT",
        help="Minimum filtering identities",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug-level logging (per-sample record counts, etc.).",
    )
    return parser.parse_args(argv)


def setup_logging(verbose: bool) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )
    return logging.getLogger("merge_abritamr")

def find_amrfinder_files(input_dir: Path, filename: str, logger: logging.Logger) -> list[Path]:
    if not input_dir.is_dir():
        logger.error(f"Input directory not found (or not a directory): {input_dir}")
        sys.exit(1)
    files = []
    for root, _dirs, filenames in os.walk(input_dir, followlinks=True):
        if filename in filenames:
            files.append(Path(root) / filename)
    files.sort()

    if not files:
        logger.error(f"No '{filename}' files found under {input_dir}")
        sys.exit(1)

    logger.info(f"Found {len(files)} '{filename}' file(s) under {input_dir}")
    return files

def load_sample(path: Path, logger: logging.Logger) -> pd.DataFrame | None:
    sample_id = path.parent.name

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", category=pd.errors.ParserWarning)
            df = pd.read_csv(
                path,
                sep="\t",
                dtype=str,
                keep_default_na=False,
                na_filter=False,
                index_col=False,
            )
        if caught:
            logger.error(
                f"[{sample_id}] {path} has a malformed row (more fields than "
                f"the header) - treating as a failed read: {caught[0].message}"
            )
            return None
    except pd.errors.EmptyDataError:
        logger.warning(f"[{sample_id}] {path} is empty (no header) - treating as a failed read")
        return None
    except pd.errors.ParserError as exc:
        logger.error(f"[{sample_id}] Could not parse {path} as TSV: {exc}")
        return None
    except OSError as exc:
        logger.error(f"[{sample_id}] Could not read {path}: {exc}")
        return None

    if df.shape[1] == 0:
        logger.warning(f"[{sample_id}] {path} has no columns - treating as a failed read")
        return None

    df.insert(0, "sample_id", sample_id)
    df.attrs["sample_id"] = sample_id

    if df.shape[0] == 0:
        logger.info(f"[{sample_id}] header-only file - 0 records")
    else:
        logger.debug(f"[{sample_id}] {df.shape[0]} record(s)")

    return df

def merge_frames(frames: list[pd.DataFrame], logger: logging.Logger) -> pd.DataFrame:
    signatures = Counter(tuple(df.columns) for df in frames)
    reference_cols, _ = signatures.most_common(1)[0]

    seen: list[str] = []
    for df in frames:
        for col in df.columns:
            if col not in seen:
                seen.append(col)
    ordered_cols = ["sample_id"] + [c for c in seen if c != "sample_id"]

    for df in frames:
        cols = tuple(df.columns)
        if cols != reference_cols:
            sample_id = df.attrs.get("sample_id", "?")
            missing = sorted(set(reference_cols) - set(cols))
            extra = sorted(set(cols) - set(reference_cols))
            logger.warning(
                f"[{sample_id}] column layout differs from the pipeline's most "
                f"common schema - missing: {missing or 'none'}, extra: {extra or 'none'}"
            )

    merged = pd.concat(frames, ignore_index=True, sort=False)
    return merged.reindex(columns=ordered_cols)
def filter_amr(
    merged: pd.DataFrame,
    output_path: Path,
    mincov_genes: float,
    minid_genes: float,
    logger: logging.Logger,
) -> None:
    TYPE_COL_IDX = 9
    COV_COL_IDX = 16
    ID_COL_IDX = 17
    required_cols = ID_COL_IDX + 1

    if merged.shape[1] < required_cols:
        logger.warning(
            f"Merged file has only {merged.shape[1]} columns (need at least "
            f"{required_cols} for the AMR/mincov/minid filter) - skipping "
            f"AMR-filtered output."
        )
        return

    type_col = merged.columns[TYPE_COL_IDX]
    cov_col = merged.columns[COV_COL_IDX]
    id_col = merged.columns[ID_COL_IDX]

    logger.info(
        f"AMR filter columns (0-indexed): [{TYPE_COL_IDX}] '{type_col}', "
        f"[{COV_COL_IDX}] '{cov_col}', [{ID_COL_IDX}] '{id_col}' "
        f"- verify these names look right for your data"
    )

    amr_filtered = merged[merged[type_col] == "AMR"].copy()

    cov_numeric = pd.to_numeric(amr_filtered[cov_col], errors="coerce")
    id_numeric = pd.to_numeric(amr_filtered[id_col], errors="coerce")

    n_bad_cov = int(cov_numeric.isna().sum())
    n_bad_id = int(id_numeric.isna().sum())
    if n_bad_cov:
        logger.warning(
            f"{n_bad_cov} AMR record(s) had a non-numeric/missing '{cov_col}' "
            f"value and were excluded by the mincov_genes filter"
        )
    if n_bad_id:
        logger.warning(
            f"{n_bad_id} AMR record(s) had a non-numeric/missing '{id_col}' "
            f"value and were excluded by the minid_genes filter"
        )

    pass_cov = cov_numeric >= mincov_genes
    pass_id = id_numeric >= minid_genes
    amr_filtered = amr_filtered[pass_cov & pass_id]

    amr_output = output_path.parent / "merged_abritamr_amr.tsv"
    amr_filtered.to_csv(amr_output, sep="\t", index=False)

    logger.info(f"AMR-filtered records      : {len(amr_filtered)} "
                f"(column '{type_col}' == 'AMR', "
                f"'{cov_col}' >= {mincov_genes}, "
                f"'{id_col}' >= {minid_genes})")
    logger.info(f"AMR output written to     : {amr_output.resolve()}")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    logger = setup_logging(args.verbose)

    files = find_amrfinder_files(args.input_dir, args.filename, logger)

    frames: list[pd.DataFrame] = []
    n_ok = 0
    n_failed = 0
    n_empty = 0

    for path in files:
        df = load_sample(path, logger)
        if df is None:
            n_failed += 1
            continue
        n_ok += 1
        if df.shape[0] == 0:
            n_empty += 1
        frames.append(df)

    if not frames:
        logger.error("No amrfinder.out files could be read successfully - nothing to merge.")
        sys.exit(1)

    merged = merge_frames(frames, logger)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    filter_amr(merged, args.output, args.mincov_genes, args.minid_genes, logger)

    logger.info("-" * 60)
    logger.info(f"amrfinder.out files found : {len(files)}")
    logger.info(f"Samples processed         : {n_ok} (failed to read: {n_failed})")
    logger.info(f"Samples with 0 records    : {n_empty}")
    logger.info(f"Total records merged      : {len(merged)}")
    logger.info(f"Output written to         : {args.output.resolve()}")
    logger.info("-" * 60)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as exc:
        logging.getLogger("merge_abritamr").error(f"Unexpected error: {exc}")
        sys.exit(1)
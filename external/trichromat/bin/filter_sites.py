#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-11 01:19


import numpy as np
import polars as pl
from scipy.stats import binom


def binom_sf_vectorized(seq):
    """Vectorized upper-tail binomial p-value P(X >= successes | trials, p).

    Equivalent to the per-row scipy.stats.binomtest(..., alternative="greater")
    but fully vectorized over all rows at once (~700x faster). Because
    scipy>=1.7 computes binomtest's "greater" p-value from the survival function,
    binom.sf(successes - 1, trials, p) is bit-identical to it.
    """
    successes = seq[0].to_numpy()
    trials = seq[1].to_numpy()
    p = seq[2].to_numpy()
    return np.where(
        (successes == 0) | (trials == 0),
        1.0,
        binom.sf(successes - 1, trials, p),
    )


def parse_and_filter(filepath, min_uncon=1, min_depth=3, min_ratio=0.05, min_pval=1):
    df_input = pl.scan_ipc(filepath)
    names = [
        c.removeprefix("Depth_")
        for c in df_input.collect_schema()
        if c.startswith("Depth_")
    ]

    df_filtered = (
        df_input.with_columns(
            (pl.col(f"Uncon_{name}") / pl.col(f"Depth_{name}")).alias(f"Ratio_{name}")
            for name in names
        )
        .with_columns(
            pl.col(f"Ratio_{name}").drop_nans().mean().alias(f"Bg_{name}")
            for name in names
        )
        .filter(pl.max_horizontal(pl.col("^Uncon_.*$")) >= min_uncon)
        .filter(pl.max_horizontal(pl.col("^Depth_.*$")) >= min_depth)
        .filter(pl.max_horizontal(pl.col("^Ratio_.*$")) >= min_ratio)
        .with_columns(
            pl.map_batches(
                [
                    pl.col(f"Uncon_{name}"),
                    pl.col(f"Depth_{name}"),
                    pl.col(f"Bg_{name}"),
                ],
                binom_sf_vectorized,
                return_dtype=pl.Float64,
            ).alias(f"pval_{name}")
            for name in names
        )
        .filter(pl.min_horizontal(pl.col("^pval_.*$")) <= min_pval)
    )
    return df_filtered


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-f", "--file", type=str, required=True, help="Path to the input file"
    )
    ap.add_argument(
        "-o", "--output", type=str, required=True, help="Path to the output file"
    )
    ap.add_argument(
        "-u",
        "--min_uncon",
        type=int,
        default=1,
        help="Minimum number of unconstrained reads",
    )
    ap.add_argument(
        "-d", "--min_depth", type=int, default=3, help="Minimum depth of coverage"
    )
    ap.add_argument(
        "-r",
        "--min_ratio",
        type=float,
        default=0.05,
        help="Minimum ratio of unconstrained reads",
    )
    ap.add_argument("-p", "--min_pval", type=float, default=1, help="Minimum p-value")

    args = ap.parse_args()

    df = parse_and_filter(
        args.file,
        min_uncon=args.min_uncon,
        min_depth=args.min_depth,
        min_ratio=args.min_ratio,
        min_pval=args.min_pval,
    )
    # polars does not support wrtting compressed file directly
    if args.output.endswith(".gz"):
        outfile = args.output[:-3]
    df.collect().write_csv(
        outfile, separator="\t", float_scientific=True, float_precision=6
    )
    if args.output.endswith(".gz"):
        import gzip
        import os
        import shutil

        with open(outfile, "rb") as f_in:
            with gzip.open(args.output, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(outfile)

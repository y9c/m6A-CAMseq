#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-10 23:33

from itertools import chain

import polars as pl


def read_file_and_rename(file, name, strand):
    df = pl.scan_ipc(file).select(
        "chrom",
        "pos",
        pl.lit(strand).alias("strand"),
        "motif",
        pl.lit(name).alias("name"),
        pl.col("unconverted").alias(f"Uncon"),
        (pl.col("converted") + pl.col("unconverted")).alias(f"Depth"),
    )
    return df


def join_table_by_site(files, names, strands):
    dfs = []
    for file, name, strand in zip(files, names, strands):
        dfs.append(read_file_and_rename(file, name, strand))

    df = (
        pl.concat(dfs)
        .group_by("chrom", "pos", "strand", "motif", maintain_order=True)
        .agg(
            pl.col(group).filter(pl.col("name") == name).sum().alias(f"{group}_{name}")
            for name in dict.fromkeys(names)
            for group in ["Uncon", "Depth"]
        )
    )
    return df


if __name__ == "__main__":
    import argparse

    arg = argparse.ArgumentParser()
    arg.add_argument(
        "-f", "--files", type=str, nargs="+", help="input files", required=True
    )
    arg.add_argument(
        "-n", "--names", type=str, nargs="+", help="input names", required=True
    )
    arg.add_argument(
        "-s", "--strands", type=str, nargs="+", help="input strands", required=True
    )
    arg.add_argument("-o", "--output", type=str, help="output file", required=True)
    args = arg.parse_args()
    # check if files and names are same length
    if len(args.files) != len(args.names):
        raise ValueError("files and names must be same length")

    df = join_table_by_site(args.files, args.names, args.strands)
    df.sink_ipc(args.output, compression="lz4")

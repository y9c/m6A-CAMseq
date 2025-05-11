#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-10 23:33

from itertools import chain

import polars as pl


def read_file_and_rename(file, name):
    df = pl.scan_csv(
        file,
        separator="\t",
        schema_overrides={
            "ref": pl.Utf8,
            "pos": pl.Int64,
            "strand": pl.Utf8,
            "convertedBaseCount": pl.Int64,
            "unconvertedBaseCount": pl.Int64,
        },
    ).select(
        "ref",
        "pos",
        "strand",
        pl.lit(name).alias("name"),
        pl.col("unconvertedBaseCount").alias(f"Uncon"),
        (pl.col("convertedBaseCount") + pl.col("unconvertedBaseCount")).alias(f"Depth"),
    )
    return df


def join_table_by_site(files, names):
    dfs = []
    for file, name in zip(files, names):
        dfs.append(read_file_and_rename(file, name))

    df = (
        pl.concat(dfs)
        .group_by("ref", "pos", "strand", maintain_order=True)
        .agg(
            pl.col(group).filter(pl.col("name") == name).sum().alias(f"{group}_{name}")
            for name in names
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
    arg.add_argument("-o", "--output", type=str, help="output file", required=True)
    args = arg.parse_args()
    # check if files and names are same length
    if len(args.files) != len(args.names):
        raise ValueError("files and names must be same length")

    df = join_table_by_site(args.files, args.names)
    df.sink_ipc(args.output, compression="lz4")

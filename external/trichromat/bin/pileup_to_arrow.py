#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-14 12:44

import polars as pl


def main(input_file: str, output_file: str):
    # Uint32 max is 4,294,967,295
    pl.scan_csv(
        input_file,
        separator="\t",
        raise_if_empty=False,
        truncate_ragged_lines=True,
        has_header=False,
        schema={
            "chrom": pl.Utf8,
            "pos": pl.UInt32,
            "motif": pl.Utf8,
            "converted": pl.UInt32,
            "unconverted": pl.UInt32,
        },
    ).sink_ipc(
        output_file,
        compression="lz4",
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert pileup to Arrow format.")
    parser.add_argument("input_file", type=str, help="Input pileup file.")
    parser.add_argument("output_file", type=str, help="Output Arrow file.")
    args = parser.parse_args()
    main(args.input_file, args.output_file)

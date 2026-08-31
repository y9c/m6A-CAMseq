#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-14 15:00

import sys
import os

def revcomp(seq):
    # str.translate is much faster for complementation
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]

class FastaIndex:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.index = {}
        self._seq = {}
        self._fh = None
        fai_file = fasta_file + ".fai"
        if not os.path.exists(fai_file):
            raise FileNotFoundError(f"FASTA index {fai_file} not found.")

        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                name = fields[0]
                length = int(fields[1])
                offset = int(fields[2])
                linebases = int(fields[3])
                linewidth = int(fields[4])
                self.index[name] = (length, offset, linebases, linewidth)

    def _load(self, chrom):
        # Read the whole chromosome once (single seek+read) instead of a
        # per-site seek+read, then cache it in memory.
        if self._fh is None:
            self._fh = open(self.fasta_file, 'r')
        length, offset, linebases, linewidth = self.index[chrom]
        nbytes = (length // linebases) * linewidth + (length % linebases)
        self._fh.seek(offset)
        raw = self._fh.read(nbytes)
        self._seq[chrom] = raw.replace("\n", "").replace("\r", "").upper()

    def fetch(self, chrom, start, end):
        if chrom not in self.index:
            return "N" * (end - start + 1)
        if chrom not in self._seq:
            self._load(chrom)

        length = self.index[chrom][0]
        seq = self._seq[chrom]

        # Clip to chromosome boundaries
        fetch_start = max(1, start)
        fetch_end = min(length, end)

        if fetch_start > fetch_end:
            return "N" * (end - start + 1)

        motif = seq[fetch_start - 1:fetch_end]

        # Pad with Ns if we clipped
        if start < 1:
            motif = "N" * (1 - start) + motif
        if end > length:
            motif = motif + "N" * (end - length)

        return motif

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Annotate motifs for hisat3n-table output.")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file.")
    parser.add_argument("-k", "--motif-flanking", type=int, default=2, help="Number of flanking bases.")
    args = parser.parse_args()

    fa = FastaIndex(args.ref)
    k = args.motif_flanking

    for line in sys.stdin:
        if line.startswith("#") or not line.strip():
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        
        chrom = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        converted = fields[3]
        unconverted = fields[4]
        
        # hisat3n-table is 1-indexed
        motif = fa.fetch(chrom, pos - k, pos + k)
        if strand == "-":
            motif = revcomp(motif)
        
        print(f"{chrom}\t{pos}\t{motif}\t{converted}\t{unconverted}")

if __name__ == "__main__":
    main()

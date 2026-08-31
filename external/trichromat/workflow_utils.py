#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-10 03:18


import os
from collections import defaultdict


def resolve_path(path, workdir):
    if os.path.isabs(path):
        return path
    else:
        path = os.path.expanduser(path)
        if os.path.isabs(path):
            return path
        else:
            return os.path.relpath(path, workdir)


def _preprocess_reference(config, workdir):
    ref = defaultdict(list)
    for k, v in config["reference"].items():
        ref[k] = [resolve_path(v2, workdir) for v2 in v]
    if "genome" not in ref:
        ref["genome"] = []
    else:
        # check if genome records is more than 1, if so, raise error
        if len(ref["genome"]) > 1:
            raise ValueError(
                "Multiple genome records found in reference. Please provide only one genome record.\n"
                "Automatically join multiple genome records into one require rebuild the index file, which takes a lot of time.\n"
            )
    # if contamination in ref but is empty, remove it
    if "contamination" in ref and len(ref["contamination"]) == 0:
        del ref["contamination"]
    return ref


def _preprocess_reads(config, workdir):
    reads = defaultdict(lambda: defaultdict(dict))
    sample2lib = defaultdict(str)
    sample2adp = defaultdict(str)
    sample2umi = defaultdict(bool)
    sample2std = defaultdict(bool)
    sample2dup = defaultdict(bool)
    sample2norc = defaultdict(bool)

    def is_umi_lib(libtype):
        return libtype in [
            "INLINE",
            "ECLIP6",
            "ECLIP10",
            "TAKARAV3",
            "SACSEQ",
            "SACSEQV3",
        ]

    def is_stranded_lib(libtype):
        return libtype not in ["UNSTRANDED"]

    for s, v in config["samples"].items():
        s = str(s)
        lib = v.get("libtype", config.get("libtype", ""))
        adp = v.get("adapter", config.get("adapter", ""))
        sample2lib[s] = lib
        sample2adp[s] = adp

        # precedence: sample > global (if exists) > inferred from libtype > default
        sample2umi[s] = v.get(
            "with_umi", config.get("with_umi", is_umi_lib(lib) if lib else True)
        )
        sample2std[s] = v.get(
            "strandness", config.get("strandness", is_stranded_lib(lib) if lib else True)
        )
        sample2dup[s] = v.get("markdup", config.get("markdup", True))
        sample2norc[s] = v.get("gene_norc", config.get("gene_norc", True))

        for i, v2 in enumerate(v["data"], 1):
            reads[s][f"run{i}"] = {
                k: resolve_path(v3, workdir) for k, v3 in dict(v2).items()
            }
    return (
        reads,
        sample2lib,
        sample2adp,
        sample2umi,
        sample2std,
        sample2dup,
        sample2norc,
    )


def preprocess_config(config, workdir):
    if "_REF" not in config:
        config["_REF"] = _preprocess_reference(config, workdir)
    if "_READS" not in config:
        (
            config["_READS"],
            config["_LIB"],
            config["_ADP"],
            config["_UMI"],
            config["_STD"],
            config["_DUP"],
            config["_NORC"],
        ) = _preprocess_reads(config, workdir)

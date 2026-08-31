#!/usr/bin/env python3

import argparse
import json
import os
import re
from collections import defaultdict

# --- Regex Patterns (for parsing file content) ---
TRIM_INPUT_READS_REGEX = re.compile(r'"input"\s*:\s*(\d+)')
TRIM_OUTPUT_READS_REGEX = re.compile(r'"output"\s*:\s*(\d+)')
READ_COUNTS_BLOCK_REGEX = re.compile(r'"read_counts"\s*:\s*{([^}]*)}')

# --- Suffixes/Patterns to strip for sample name extraction, by file type context ---
SAMPLE_NAME_STRIP_PATTERNS = {
    "trimming_json": re.compile(r"_run\d+\.json$"),
    "combined_genome_count": re.compile(r"\.genome\.count$"),
    "combined_genes_count": re.compile(r"\.genes\.count$"),
    "combined_contamination_count": re.compile(r"\.contamination\.count$"),
    "dedup_genome_count": re.compile(r"\.genome\.count$"),
    "dedup_genes_count": re.compile(r"\.genes\.count$"),
    "dedup_contamination_count": re.compile(r"\.contamination\.count$"),
    "unmap_count": re.compile(r"\.count$"),
}

# Regex to determine type (genome, genes, contamination) from .count filenames
FILENAME_TYPE_REGEX = re.compile(r"\.(genome|genes|contamination)\.count$")


def get_sample_name(file_path, file_type_context):
    """
    Extracts a sample name from a file path by stripping known suffixes
    based on file_type_context.
    """
    base_name = os.path.basename(file_path)

    pattern_to_strip = SAMPLE_NAME_STRIP_PATTERNS.get(file_type_context)

    if pattern_to_strip:
        sample_name = pattern_to_strip.sub("", base_name)
        if sample_name != base_name and sample_name:
            return sample_name

    fallback_name = base_name.split(".", 1)[0]
    if not fallback_name:
        fallback_name = base_name
    print(
        f"Warning: Using fallback sample name extraction for '{base_name}' (context: '{file_type_context}'). Extracted: '{fallback_name}'"
    )
    return fallback_name


def parse_trimming_json(file_path):
    """Parses a trimming JSON report for input and output read counts."""
    try:
        with open(file_path, "r") as f:
            content = f.read()
        try:
            data = json.loads(content)
            rc = data.get("read_counts", {})
            inp, outp = rc.get("input"), rc.get("output")
            if inp is not None and outp is not None:
                return {"raw": int(inp), "clean": int(outp)}
        except json.JSONDecodeError:
            print(f"Warning: JSON decode error for {file_path}. Falling back to regex.")

        rb_match = READ_COUNTS_BLOCK_REGEX.search(content)
        if rb_match:
            rc_content = rb_match.group(1)
            inp_match = TRIM_INPUT_READS_REGEX.search(rc_content)
            outp_match = TRIM_OUTPUT_READS_REGEX.search(rc_content)
            if inp_match and outp_match:
                return {
                    "raw": int(inp_match.group(1)),
                    "clean": int(outp_match.group(1)),
                }
        print(f"Warning: Could not parse trimming stats from {file_path}.")
        return {"raw": None, "clean": None}
    except Exception as e:
        print(f"Error parsing trimming file {file_path}: {e}")
        return {"raw": None, "clean": None}


def parse_simple_count_file(file_path):
    """Parses a file assumed to contain a single number."""
    try:
        with open(file_path, "r") as f:
            content = f.read().strip()
        if content:
            return int(content)
        print(f"Warning: Simple count file {file_path} is empty.")
        return None
    except ValueError:
        print(f"Warning: Could not convert content of {file_path} to int: '{content}'")
        return None
    except Exception as e:
        print(f"Error parsing simple count file {file_path}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate read counts from various bioinformatics report files.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--output-name",
        default="read_counts_summary.tsv",
        help="Name for the output TSV file. If extension is not .tsv, it will be appended. Default: 'read_counts_summary.tsv'",
    )

    parser.add_argument(
        "--trimming", nargs="+", help="Trimming JSON reports (e.g., Sample_runX.json)."
    )

    parser.add_argument(
        "--unmap",
        nargs="+",
        required=False,
        help="Unmapped read count files (e.g., Sample.count).",
    )

    parser.add_argument(
        "--combined",
        nargs="+",
        help="Mapping results from 'combined' dir (e.g., Sample.genome.count, Sample.genes.count).",
    )
    parser.add_argument(
        "--dedup",
        nargs="+",
        help="Final deduplicated counts from 'dedup' dir (e.g., Sample.genome.count, Sample.genes.count).",
    )

    args = parser.parse_args()
    results = defaultdict(lambda: defaultdict(lambda: None))

    # --- Process Trimming ---
    if args.trimming:
        for fp in args.trimming:
            s_name = get_sample_name(fp, "trimming_json")
            counts = parse_trimming_json(fp)
            if counts["raw"] is not None:
                results[s_name]["Raw_Reads"] = (
                    results[s_name].get("Raw_Reads", 0) or 0
                ) + counts["raw"]
            if counts["clean"] is not None:
                results[s_name]["Clean_Reads"] = (
                    results[s_name].get("Clean_Reads", 0) or 0
                ) + counts["clean"]

    # --- Process Unmap ---
    if args.unmap:
        for fp in args.unmap:
            s_name = get_sample_name(fp, "unmap_count")
            results[s_name]["Unmap_Reads"] = parse_simple_count_file(fp)

    # --- Process Combined Directory Results ---
    if args.combined:
        for fp in args.combined:
            base_name = os.path.basename(fp)
            type_match = FILENAME_TYPE_REGEX.search(base_name)
            if type_match and fp.endswith(".count"):
                file_type = type_match.group(1)  # genome, genes, or contamination
                s_name_context = f"combined_{file_type}_count"
                s_name = get_sample_name(fp, s_name_context)
                metric_name = f"Mapped_{file_type.capitalize()}_Reads"
                results[s_name][metric_name] = parse_simple_count_file(fp)
            else:
                print(
                    f"Warning: File '{fp}' for --combined could not be categorized or is not a .count file. Expected format like 'Sample.type.count'. Skipping."
                )

    # --- Process Dedup Results ---
    if args.dedup:
        for fp in args.dedup:
            base_name = os.path.basename(fp)
            type_match = FILENAME_TYPE_REGEX.search(base_name)
            if type_match and fp.endswith(".count"):
                file_type = type_match.group(1)  # genome, genes, or contamination
                s_name_context = f"dedup_{file_type}_count"
                s_name = get_sample_name(fp, s_name_context)
                metric_name = f"Dedup_{file_type.capitalize()}_Reads"
                results[s_name][metric_name] = parse_simple_count_file(fp)
            else:
                print(
                    f"Warning: File '{fp}' for --dedup could not be categorized or is not a .count file. Expected format like 'Sample.type.count'. Skipping."
                )

    # --- Output TSV Matrix ---
    if not results:
        print("No data processed. Exiting.")
        return

    output_filename = args.output_name
    if not output_filename.lower().endswith(".tsv"):
        output_filename += ".tsv"

    header = ["Sample"]
    all_metrics = set()
    for data in results.values():
        all_metrics.update(data.keys())

    final_preferred_order = [
        "Raw_Reads",
        "Clean_Reads",
        "Unmap_Reads",
        "Mapped_Contamination_Reads",
        "Mapped_Genes_Reads",
        "Mapped_Genome_Reads",
        "Dedup_Contamination_Reads",
        "Dedup_Genes_Reads",
        "Dedup_Genome_Reads",
    ]

    current_header_metrics = [m for m in final_preferred_order if m in all_metrics]
    current_header_metrics.extend(
        sorted([m for m in all_metrics if m not in current_header_metrics])
    )
    header.extend(current_header_metrics)

    print(f"\nWriting aggregated results to: {output_filename}")
    try:
        with open(output_filename, "w") as outfile:
            outfile.write("\t".join(header) + "\n")
            for s_name in sorted(results.keys()):
                row_values = [s_name]
                for metric_name in header[1:]:
                    value = results[s_name].get(metric_name)
                    row_values.append(str(value) if value is not None else "NA")
                outfile.write("\t".join(row_values) + "\n")
        print(f"Successfully wrote summary to {output_filename}")
    except Exception as e:
        print(f"Error writing output file {output_filename}: {e}")
    print("Script finished.")


if __name__ == "__main__":
    main()

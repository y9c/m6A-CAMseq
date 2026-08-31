#include <ctype.h>  // For toupper, isdigit, isalpha
#include <errno.h>  // For EPIPE
#include <getopt.h> // For command-line option parsing
#include <htslib/faidx.h>
#include <htslib/hts.h>     // For hts_pos_t and other common htslib types
#include <htslib/kstring.h> // For kstring_t
#include <htslib/sam.h>
#include <inttypes.h> // For PRId64 if used in logging
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // For PG line date

// Define conversion/unconversion characters based on strand
typedef struct {
  char r0_ref;  // Expected reference base for primary conversion
  char a0_conv; // Observed read base if primary conversion occurred
  char r1_ref;  // Expected reference base for secondary conversion
  char a1_conv; // Observed read base if secondary conversion occurred
} conversion_defs_t;

// Global variables for user-defined mutations (defaults for '+' strand)
char primary_ref_char_plus_strand = 'A';
char primary_alt_char_plus_strand = 'G';
char secondary_ref_char_plus_strand = 'C';
char secondary_alt_char_plus_strand = 'T';

// Helper to get complement of a base
char complement_base(char base) {
  base = toupper(base);
  switch (base) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'N':
    return 'N';
  default:
    return 'N';
  }
}

// Helper to convert 4-bit sequence to char
#define base_to_char(b) (seq_nt16_str[(b)])

// Function to determine effective strand and conversion definitions
char determine_strand_and_defs(const bam1_t *aln, int reversed_lib_flag,
                               conversion_defs_t *defs) {
  const bam1_core_t *c = &aln->core;
  char effective_strand_char;

  if (c->flag & BAM_FPAIRED) {
    int is_read1 = (c->flag & BAM_FREAD1) != 0;
    int is_reverse_mapped = (c->flag & BAM_FREVERSE) != 0;
    if ((is_read1 && is_reverse_mapped) || (!is_read1 && !is_reverse_mapped)) {
      effective_strand_char = '-';
    } else {
      effective_strand_char = '+';
    }
  } else {
    effective_strand_char = (c->flag & BAM_FREVERSE) ? '-' : '+';
  }

  if (reversed_lib_flag) {
    effective_strand_char = (effective_strand_char == '+') ? '-' : '+';
  }

  if (effective_strand_char == '+') {
    defs->r0_ref = primary_ref_char_plus_strand;
    defs->a0_conv = primary_alt_char_plus_strand;
    defs->r1_ref = secondary_ref_char_plus_strand;
    defs->a1_conv = secondary_alt_char_plus_strand;
  } else { // '-' strand
    defs->r0_ref = complement_base(primary_ref_char_plus_strand);
    defs->a0_conv = complement_base(primary_alt_char_plus_strand);
    defs->r1_ref = complement_base(secondary_ref_char_plus_strand);
    defs->a1_conv = complement_base(secondary_alt_char_plus_strand);
  }
  return effective_strand_char;
}

// Process a single BAM record and add tags
int process_and_tag_read(bam1_t *aln, bam_hdr_t *header, faidx_t *fai,
                         int reversed_lib_flag) {
  bam1_core_t *c = &aln->core;

  int yf = 0, zf = 0, yc = 0, zc = 0, ns = 0, nc = 0;
  char *ref_segment_fai = NULL;
  hts_pos_t ref_segment_len_fai = 0;

  if (c->flag & BAM_FUNMAP) {
    nc = c->l_qseq > 0 ? c->l_qseq : 0;
    goto set_tags_and_return_unmapped;
  }

  conversion_defs_t defs;
  determine_strand_and_defs(aln, reversed_lib_flag, &defs);

  uint8_t *query_seq = bam_get_seq(aln);
  uint32_t *cigar = bam_get_cigar(aln);
  const char *md_z = NULL;
  uint8_t *md_aux = bam_aux_get(aln, "MD");

  if (md_aux && (*md_aux == 'Z' || *md_aux == 'H')) {
    md_z = bam_aux2Z(md_aux);
  }

  hts_pos_t ref_fetch_start = c->pos;
  hts_pos_t ref_fetch_end = bam_endpos(aln) - 1;

  if (fai && c->tid >= 0 && header->target_name &&
      header->target_name[c->tid] && ref_fetch_end >= ref_fetch_start) {
    ref_segment_fai =
        faidx_fetch_seq64(fai, header->target_name[c->tid], ref_fetch_start,
                          ref_fetch_end, &ref_segment_len_fai);
    if (!ref_segment_fai) {
      // fprintf(stderr, "Warning: Could not fetch FASTA for read %s
      // (%s:%"PRId64"-%"PRId64") L=%lld. MD tag is crucial.\n",
      //         bam_get_qname(aln), header->target_name[c->tid],
      //         ref_fetch_start, ref_fetch_end, (long
      //         long)ref_segment_len_fai);
    }
  }

  if (!md_z && !ref_segment_fai) {
    nc = c->l_qseq > 0 ? c->l_qseq : 0;
    goto set_tags_and_return_info_missing;
  }

  int q_idx = 0;
  hts_pos_t r_idx_fai_offset = 0;

  const char *md_ptr = md_z;
  int md_match_countdown = 0;

  for (uint32_t i = 0; i < c->n_cigar; ++i) {
    int cigar_op = bam_cigar_op(cigar[i]);
    int cigar_len = bam_cigar_oplen(cigar[i]);

    if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL ||
        cigar_op == BAM_CDIFF) {
      for (int k = 0; k < cigar_len; ++k) {
        char read_base_char = 'N';
        if (q_idx < c->l_qseq) {
          read_base_char = toupper(base_to_char(bam_seqi(query_seq, q_idx)));
        } else {
          nc++;
          q_idx++;
          r_idx_fai_offset++;
          continue;
        }

        char ref_base_char = 'N';
        if (md_ptr) {
          if (md_match_countdown > 0) {
            ref_base_char = read_base_char;
            md_match_countdown--;
          } else {
            if (isdigit((unsigned char)*md_ptr)) {
              char *end_num_ptr;
              md_match_countdown = strtol(md_ptr, &end_num_ptr, 10);
              md_ptr = end_num_ptr;
              if (md_match_countdown > 0) {
                ref_base_char = read_base_char;
                md_match_countdown--;
              }
            }
            if (ref_base_char == 'N' && isalpha((unsigned char)*md_ptr)) {
              ref_base_char = toupper((unsigned char)*md_ptr);
              md_ptr++;
            } else if (ref_base_char == 'N' && *md_ptr == '^') {
              ref_base_char = 'N';
            }
          }
        } else if (ref_segment_fai && r_idx_fai_offset < ref_segment_len_fai) {
          ref_base_char = toupper(ref_segment_fai[r_idx_fai_offset]);
        }

        // ** CORRECTED LOGIC FOR ns **
        if (read_base_char == 'N' || ref_base_char == 'N') {
          nc++;
        } else if (ref_base_char == defs.r0_ref) { // Primary site
          if (read_base_char == defs.a0_conv) {    // Primary conversion
            yf++;
          } else if (read_base_char ==
                     defs.r0_ref) { // Primary unconverted (match)
            zf++;
          } else { // Substitution at a primary site
            ns++;
          }
        } else if (ref_base_char == defs.r1_ref) { // Secondary site
          if (read_base_char == defs.a1_conv) {    // Secondary conversion
            yc++;
          } else if (read_base_char ==
                     defs.r1_ref) { // Secondary unconverted (match)
            zc++;
          } else { // Substitution at a secondary site
            ns++;
          }
        } else { // Neither primary nor secondary reference base
          if (read_base_char != ref_base_char) { // Mismatch at an "other" site
            ns++;
          }
          // If read_base_char == ref_base_char here, it's a match at an "other"
          // site, no specific counter (yf,zf,yc,zc,ns) incremented.
        }
        q_idx++;
        r_idx_fai_offset++;
      }
    } else if (cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP) {
      for (int k = 0; k < cigar_len; ++k) {
        nc++;
        q_idx++;
      }
    } else if (cigar_op == BAM_CDEL) {
      if (md_ptr && *md_ptr == '^') {
        md_ptr++;
        for (int k = 0; k < cigar_len; ++k) {
          nc++;
          if (isalpha((unsigned char)*md_ptr))
            md_ptr++;
          else
            break;
          r_idx_fai_offset++;
        }
      } else {
        for (int k = 0; k < cigar_len; ++k) {
          nc++;
          r_idx_fai_offset++;
        }
      }
    } else if (cigar_op == BAM_CREF_SKIP || cigar_op == BAM_CPAD) {
      for (int k = 0; k < cigar_len; ++k) {
        nc++;
        r_idx_fai_offset++;
      }
    }
  }

  // Common exit point for successfully processed mapped reads
  bam_aux_update_int(aln, "Yf", yf);
  bam_aux_update_int(aln, "Zf", zf);
  bam_aux_update_int(aln, "Yc", yc);
  bam_aux_update_int(aln, "Zc", zc);
  bam_aux_update_int(aln, "NS", ns);
  bam_aux_update_int(aln, "NC", nc);
  if (ref_segment_fai)
    free(ref_segment_fai);
  return 0;

set_tags_and_return_unmapped:
  bam_aux_update_int(aln, "Yf", 0);
  bam_aux_update_int(aln, "Zf", 0);
  bam_aux_update_int(aln, "Yc", 0);
  bam_aux_update_int(aln, "Zc", 0);
  bam_aux_update_int(aln, "NS", 0);
  bam_aux_update_int(aln, "NC", nc);
  // ref_segment_fai is NULL here, no need to free
  return 1;

set_tags_and_return_info_missing:
  bam_aux_update_int(aln, "Yf", 0);
  bam_aux_update_int(aln, "Zf", 0);
  bam_aux_update_int(aln, "Yc", 0);
  bam_aux_update_int(aln, "Zc", 0);
  bam_aux_update_int(aln, "NS", 0);
  bam_aux_update_int(aln, "NC", nc);
  if (ref_segment_fai)
    free(ref_segment_fai);
  return 2;
}

void print_usage(const char *progname) {
  fprintf(stderr, "Usage: %s [options] [input.bam [output.bam]]\n\n", progname);
  fprintf(stderr, "If input.bam is '-' or missing, reads from stdin.\n");
  fprintf(stderr,
          "If output.bam is '-' or missing, writes to stdout (as BAM).\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr,
          "  -r, --reference <ref.fasta>   Reference FASTA file (indexed).\n");
  fprintf(stderr, "  -p, --primary-mutation <RF>   Define primary mutation for "
                  "'+' strand (e.g., AG for A->G). Default: AG\n");
  fprintf(stderr, "  -s, --secondary-mutation <RF> Define secondary mutation "
                  "for '+' strand (e.g., CT for C->T). Default: CT\n");
  fprintf(stderr,
          "  -L, --reversed-lib            Treat library as reversed type "
          "(strand interpretation is flipped for conversion defs).\n");
  fprintf(stderr, "                                By default, library is "
                  "treated as forward.\n");
  fprintf(stderr,
          "  -h, --help                      Show this help message.\n\n");
  fprintf(stderr, "Description: Calculates conversion counts based on "
                  "reference/query comparison and\n");
  fprintf(stderr,
          "             adds Yf, Zf, Yc, Zc, NS, NC tags to each read.\n");
  fprintf(stderr, "             Requires MD tag in input BAM or a reference "
                  "FASTA file (MD tag prioritized).\n");
}

int main(int argc, char *argv[]) {
  char *bam_file_in_path = "-";
  char *bam_file_out_path = "-";
  int reversed_lib_cli_val = 0;
  char *ref_fasta_path = NULL;
  char *primary_mutation_str = "AG";
  char *secondary_mutation_str = "CT";

  int opt;
  static struct option long_options[] = {
      {"reference", required_argument, 0, 'r'},
      {"primary-mutation", required_argument, 0, 'p'},
      {"secondary-mutation", required_argument, 0, 's'},
      {"reversed-lib", no_argument, 0, 'L'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};

  kstring_t pg_cl = {0, 0, 0};
  for (int i = 0; i < argc; ++i) {
    if (pg_cl.l > 0)
      kputc(' ', &pg_cl);
    int has_space = 0;
    for (char *p = argv[i]; *p; p++)
      if (isspace((unsigned char)*p)) {
        has_space = 1;
        break;
      } // cast to unsigned char
    if (has_space)
      kputc('"', &pg_cl);
    kputs(argv[i], &pg_cl);
    if (has_space)
      kputc('"', &pg_cl);
  }

  while ((opt = getopt_long(argc, argv, "hr:p:s:L", long_options, NULL)) !=
         -1) {
    switch (opt) {
    case 'h':
      print_usage(argv[0]);
      free(pg_cl.s);
      return 0;
    case 'r':
      ref_fasta_path = optarg;
      break;
    case 'p':
      primary_mutation_str = optarg;
      break;
    case 's':
      secondary_mutation_str = optarg;
      break;
    case 'L':
      reversed_lib_cli_val = 1;
      break;
    case '?':
      print_usage(argv[0]);
      free(pg_cl.s);
      return 1;
    default:
      abort();
    }
  }

  if (optind < argc) {
    bam_file_in_path = argv[optind++];
  }
  if (optind < argc) {
    bam_file_out_path = argv[optind++];
  }
  if (optind < argc) {
    fprintf(stderr, "Error: Too many positional arguments.\n");
    print_usage(argv[0]);
    free(pg_cl.s);
    return 1;
  }

  if (strlen(primary_mutation_str) == 2 &&
      isalpha((unsigned char)primary_mutation_str[0]) &&
      isalpha((unsigned char)primary_mutation_str[1])) {
    primary_ref_char_plus_strand = toupper(primary_mutation_str[0]);
    primary_alt_char_plus_strand = toupper(primary_mutation_str[1]);
  } else {
    fprintf(stderr,
            "Error: Primary mutation string '%s' is invalid. Must be two "
            "alphabetic characters (e.g., AG).\n",
            primary_mutation_str);
    free(pg_cl.s);
    return 1;
  }

  if (strlen(secondary_mutation_str) == 2 &&
      isalpha((unsigned char)secondary_mutation_str[0]) &&
      isalpha((unsigned char)secondary_mutation_str[1])) {
    secondary_ref_char_plus_strand = toupper(secondary_mutation_str[0]);
    secondary_alt_char_plus_strand = toupper(secondary_mutation_str[1]);
  } else {
    fprintf(stderr,
            "Error: Secondary mutation string '%s' is invalid. Must be two "
            "alphabetic characters (e.g., CT).\n",
            secondary_mutation_str);
    free(pg_cl.s);
    return 1;
  }

  // Only print to stderr if not writing primary output to stdout, to avoid
  // clutter
  if (strcmp(bam_file_out_path, "-") != 0) {
    fprintf(stderr, "INFO: Input source: %s\n",
            strcmp(bam_file_in_path, "-") == 0 ? "stdin" : bam_file_in_path);
    fprintf(stderr, "INFO: Output destination: %s\n",
            bam_file_out_path); // Already known not to be stdout
    fprintf(stderr, "INFO: Primary mutation (+ strand): %c -> %c\n",
            primary_ref_char_plus_strand, primary_alt_char_plus_strand);
    fprintf(stderr, "INFO: Secondary mutation (+ strand): %c -> %c\n",
            secondary_ref_char_plus_strand, secondary_alt_char_plus_strand);
    fprintf(stderr, "INFO: Library type: %s\n",
            reversed_lib_cli_val ? "Reversed" : "Forward (default)");
  }

  faidx_t *fai = NULL;
  unsigned long long reads_written_total = 0;
  unsigned long long reads_processed_mapped_ok = 0;
  unsigned long long reads_processed_unmapped = 0;
  unsigned long long reads_processed_info_missing = 0;

  if (ref_fasta_path) {
    fai = fai_load(ref_fasta_path);
    if (!fai) {
      fprintf(stderr,
              "Warning: Could not load FASTA index for %s. MD tag will be "
              "essential.\n",
              ref_fasta_path);
    } else {
      if (strcmp(bam_file_out_path, "-") != 0)
        fprintf(stderr, "INFO: Loaded reference FASTA: %s\n", ref_fasta_path);
    }
  } else {
    if (strcmp(bam_file_out_path, "-") != 0)
      fprintf(stderr,
              "INFO: No reference FASTA provided. MD tag will be essential.\n");
  }

  htsFile *fp_in = NULL, *fp_out = NULL;
  bam_hdr_t *header = NULL;
  bam1_t *aln = NULL;
  int ret_read, ret_write;
  const char *in_mode = "r";
  const char *out_mode = "wb";

  if (strcmp(bam_file_in_path, "-") == 0 &&
      strcmp(bam_file_out_path, "-") != 0) {
    fprintf(stderr, "INFO: Reading from stdin.\n");
  }
  if (strcmp(bam_file_out_path, "-") == 0 &&
      strcmp(bam_file_in_path, "-") != 0) {
    fprintf(stderr, "INFO: Writing BAM output to stdout.\n");
  }
  if (strcmp(bam_file_out_path, "-") == 0 &&
      strcmp(bam_file_in_path, "-") == 0) {
    fprintf(stderr,
            "INFO: Reading from stdin and writing BAM output to stdout.\n");
  }

  if ((fp_in = hts_open(bam_file_in_path, in_mode)) == NULL) {
    fprintf(stderr, "Error: Could not open input source '%s': %s\n",
            bam_file_in_path, strerror(errno));
    if (fai)
      fai_destroy(fai);
    free(pg_cl.s);
    return 1;
  }

  if ((header = sam_hdr_read(fp_in)) == NULL) {
    fprintf(stderr, "Error: Could not read header from '%s'\n",
            bam_file_in_path);
    hts_close(fp_in);
    if (fai)
      fai_destroy(fai);
    free(pg_cl.s);
    return 1;
  }

  char version_str[64];
  time_t t = time(NULL);
  struct tm tm_info;
  localtime_r(&t, &tm_info);
  snprintf(version_str, sizeof(version_str), "1.0.4-%04d%02d%02d",
           tm_info.tm_year + 1900, tm_info.tm_mon + 1, tm_info.tm_mday);

  if (sam_hdr_add_pg(header, "add_conversion_tags", "VN", version_str, "CL",
                     pg_cl.s ? pg_cl.s : "", NULL) < 0) {
    fprintf(stderr, "Warning: Failed to add PG line to header.\n");
  }
  free(pg_cl.s);

  if ((fp_out = hts_open(bam_file_out_path, out_mode)) == NULL) {
    fprintf(stderr, "Error: Could not open output destination '%s': %s\n",
            bam_file_out_path, strerror(errno));
    sam_hdr_destroy(header);
    hts_close(fp_in);
    if (fai)
      fai_destroy(fai);
    return 1;
  }

  if (sam_hdr_write(fp_out, header) < 0) {
    fprintf(stderr, "Error: Could not write header to '%s'\n",
            bam_file_out_path);
    hts_close(fp_out);
    sam_hdr_destroy(header);
    hts_close(fp_in);
    if (fai)
      fai_destroy(fai);
    return 1;
  }

  if ((aln = bam_init1()) == NULL) {
    fprintf(stderr, "Error: Could not initialize BAM alignment structure.\n");
    hts_close(fp_out);
    sam_hdr_destroy(header);
    hts_close(fp_in);
    if (fai)
      fai_destroy(fai);
    return 1;
  }

  if (strcmp(bam_file_out_path, "-") != 0)
    fprintf(stderr, "INFO: Processing reads...\n");
  while ((ret_read = sam_read1(fp_in, header, aln)) >= 0) {
    int process_status =
        process_and_tag_read(aln, header, fai, reversed_lib_cli_val);

    if (process_status == 1) {
      reads_processed_unmapped++;
    } else if (process_status == 2) {
      reads_processed_info_missing++;
    } else {
      reads_processed_mapped_ok++;
    }

    ret_write = sam_write1(fp_out, header, aln);
    if (ret_write < 0) {
      fprintf(stderr, "Error: Failed to write BAM record to '%s'.\n",
              bam_file_out_path);
      if (strcmp(bam_file_out_path, "-") == 0 &&
          (errno == EPIPE || ferror(stdout))) { // Check ferror for stdout
        fprintf(stderr, "       (Broken pipe or error writing to stdout - was "
                        "stdout closed by downstream process?)\n");
      } else if (strcmp(bam_file_out_path, "-") != 0) {
        fprintf(stderr, "       %s\n", strerror(errno));
      }
      break;
    }
    reads_written_total++;
  }
  if (ret_read < -1) {
    fprintf(stderr,
            "Warning: Error reading input BAM record from '%s' (sam_read1 "
            "returned %d).\n",
            bam_file_in_path, ret_read);
  }
  if (strcmp(bam_file_out_path, "-") != 0) {
    fprintf(stderr, "INFO: Finished processing. Total reads written: %llu\n",
            reads_written_total);

    fprintf(stderr, "--- Processing Summary ---\n");
    fprintf(stderr, "Input Source: %s\n",
            strcmp(bam_file_in_path, "-") == 0 ? "stdin" : bam_file_in_path);
    fprintf(stderr, "Output Destination: %s\n",
            strcmp(bam_file_out_path, "-") == 0 ? "stdout" : bam_file_out_path);
    fprintf(stderr, "Library Type: %s\n",
            reversed_lib_cli_val ? "Reversed" : "Forward (default)");
    if (ref_fasta_path)
      fprintf(stderr, "Reference FASTA: %s %s\n", ref_fasta_path,
              fai ? "(Loaded)" : "(Not loaded successfully)");
    else
      fprintf(stderr, "Reference FASTA: Not provided.\n");
    fprintf(stderr, "Primary Mutation (+strand): %c->%c\n",
            primary_ref_char_plus_strand, primary_alt_char_plus_strand);
    fprintf(stderr, "Secondary Mutation (+strand): %c->%c\n",
            secondary_ref_char_plus_strand, secondary_alt_char_plus_strand);

    fprintf(stderr, "Total Reads Written to Output: %llu\n",
            reads_written_total);
    fprintf(stderr, "  - Mapped & successfully processed: %llu\n",
            reads_processed_mapped_ok);
    fprintf(stderr, "  - Unmapped (default tags): %llu\n",
            reads_processed_unmapped);
    fprintf(stderr, "  - Mapped but missing MD/ref info (default tags): %llu\n",
            reads_processed_info_missing);
    fprintf(stderr,
            "--------------------------------------------------------\n");
  }

  bam_destroy1(aln);
  sam_hdr_destroy(header);
  if (fp_in && hts_close(fp_in) < 0) {
    if (strcmp(bam_file_in_path, "-") !=
        0) /* Don't report error for closing stdin */
      fprintf(stderr, "Warning: Error closing input source '%s'.\n",
              bam_file_in_path);
  }
  if (fp_out && hts_close(fp_out) < 0) {
    if (strcmp(bam_file_out_path, "-") !=
        0) /* Don't report error for closing stdout */
      fprintf(stderr, "Warning: Error closing output destination '%s'.\n",
              bam_file_out_path);
  }
  if (fai)
    fai_destroy(fai);

  return 0;
}

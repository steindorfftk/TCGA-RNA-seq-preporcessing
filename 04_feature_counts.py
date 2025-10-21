import os
import argparse
import time
from pathlib import Path
from utils import run_command

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="""
Run featureCounts on STAR-aligned BAM files.

This script counts reads mapping to genes using featureCounts.
It supports single-end (SE) and paired-end (PE) reads, and can optionally
delete BAM files after processing to save memory.
"""
    )

    parser.add_argument(
        "-nv", "--verbose", action="store_true", default=True,
        help="Enable verbose output (default: True)."
    )
    parser.add_argument(
        "-e", "--endtype", choices=["se", "pe"], required=True,
        help="Specify sequencing type: 'se' for single-end or 'pe' for paired-end reads."
    )
    parser.add_argument(
        "-lm", "--lowmemory", action="store_true",
        help="Delete BAM files after processing to save memory."
    )

    args = parser.parse_args()

    # Annotation GTF file
    annotation_path = "temporary/star/annotation/gencode.v36.annotation.gtf"

    # Directories
    star_output_dir = Path("temporary/star/second_pass_output")
    counts_output_dir = Path("output")
    os.makedirs(counts_output_dir, exist_ok=True)

    # Locate BAM files
    bam_files = sorted(star_output_dir.glob("*Aligned.sortedByCoord.out.bam"))
    if not bam_files:
        print("❌ No BAM files found in temporary/star/second_pass_output/. Check your STAR run.")
        return

    # featureCounts flag for paired-end
    paired_flag = "-p" if args.endtype == "pe" else ""

    for bam_file in bam_files:
        srr = bam_file.stem.split("_")[0]  # extract SRR ID
        output_file = counts_output_dir / f"{srr}_featurecounts.txt"

        # Construct featureCounts command
        command = (
            f"featureCounts {paired_flag} -T 8 "
            f"-t exon -g gene_id "
            f"-a {annotation_path} "
            f"-o {output_file} {bam_file}"
        )

        run_command(command, args.verbose)

        if args.verbose:
            print(f"✅ Feature counting complete for {srr}\n")

    # Delete BAM files if low-memory mode is enabled
    if args.lowmemory:
        for bam_file in bam_files:
            bam_file.unlink()
        if args.verbose:
            print("🧹 BAM files deleted to save memory.\n")

    # Clean and extract gene counts
    for txt_file in counts_output_dir.glob("*_featurecounts.txt"):
        srr = txt_file.stem.split("_")[0]
        clean_output = counts_output_dir / f"{srr}_gene_counts.tabular"
        os.makedirs(clean_output.parent, exist_ok=True)

        # Remove header lines and keep only gene_id + count column (column 7)
        run_command(f"sed -i '/^#/d' {txt_file}", args.verbose)
        run_command(f"cut -f 1,7 {txt_file} > {clean_output}", args.verbose)

    elapsed_time = time.time() - start_time
    if args.verbose:
        print(f"✅ FeatureCounts completed successfully. Total runtime: {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()


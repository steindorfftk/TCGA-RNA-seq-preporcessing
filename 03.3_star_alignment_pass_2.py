import os
import argparse
import time
from utils import run_command

def main():
    parser = argparse.ArgumentParser(
        description="""
Perform STAR second-pass RNA-seq alignment.

This script aligns reads using the second-pass STAR genome index to improve mapping accuracy.
Supports single-end (SE) and paired-end (PE) reads, and optionally uses trimmed FASTQs.

Workflow:
1. Reads accession numbers from input/SRRAccList.txt
2. Aligns reads to the second-pass genome index
3. Outputs sorted BAM files and alignment statistics
"""
    )

    parser.add_argument(
        "-o", "--organism", choices=["Homo", "Mus"], required=True,
        help="Specify organism: Homo (human) or Mus (mouse)."
    )
    parser.add_argument(
        "-e", "--ending", choices=["se", "pe"], required=True,
        help="Specify sequencing type: single-end (SE) or paired-end (PE)."
    )
    parser.add_argument(
        "-nv", "--verbose", action="store_true", default=True,
        help="Enable verbose output (default: True)."
    )
    parser.add_argument(
        "-t", "--trimmed", action="store_true",
        help="Use trimmed FASTQs from Trimmomatic instead of raw FASTQs."
    )
    parser.add_argument(
        "-th", "--threads", type=int, default=6,
        help="Number of threads to use for STAR (default: 6)."
    )

    args = parser.parse_args()

    # Set directories
    input_dir = "temporary/trimmomatic/output" if args.trimmed else "temporary/fastq_dump"
    second_pass_index_dir = "temporary/star/index/GRCh38_d1_vd1_v36_secondpass"
    second_pass_output_dir = "temporary/star/second_pass_output"
    stats_dir = os.path.join(second_pass_output_dir, "stats")

    os.makedirs(second_pass_output_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    # Read SRR accession list
    srr_list = [line.strip() for line in open("input/SRRAccList.txt")]

    print("=== STEP 4: Second-pass STAR alignment ===")
    for srr in srr_list:
        start_time = time.time()

        # Prepare input FASTQ files
        if args.ending == "pe":
            read_files_in = f"{os.path.join(input_dir, srr+'_1.fastq')} {os.path.join(input_dir, srr+'_2.fastq')}"
        else:
            read_files_in = os.path.join(input_dir, srr + ".fastq")

        output_prefix = os.path.join(second_pass_output_dir, f"{srr}_")
        rg_line = f"ID:{srr} PL:ILLUMINA LB:{args.organism} SM:{srr}"

        # Build STAR command
        command = (
            f"STAR "
            f"--genomeDir {second_pass_index_dir} "
            f"--readFilesIn {read_files_in} "
            f"--runThreadN {args.threads} "
            f"--outFileNamePrefix {output_prefix} "
            f"--outFilterMultimapScoreRange 1 "
            f"--outFilterMultimapNmax 20 "
            f"--outFilterMismatchNmax 10 "
            f"--alignIntronMax 500000 "
            f"--alignMatesGapMax 1000000 "
            f"--sjdbScore 2 "
            f"--alignSJDBoverhangMin 1 "
            f"--genomeLoad NoSharedMemory "
            f"--limitBAMsortRAM 0 "
            f"--readFilesCommand cat "
            f"--outFilterMatchNminOverLread 0.33 "
            f"--outFilterScoreMinOverLread 0.33 "
            f"--sjdbOverhang 100 "
            f"--outSAMstrandField intronMotif "
            f"--outSAMattributes NH HI NM MD AS XS "
            f"--outSAMunmapped Within "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outSAMheaderHD '@HD VN:1.4' "
            f"--outSAMattrRGline '{rg_line}' "
            f"2> {stats_dir}/{srr}_secondpass_stats.txt"
        )

        # Execute alignment
        run_command(command, args.verbose)
        print(f"Completed second-pass for {srr} in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()


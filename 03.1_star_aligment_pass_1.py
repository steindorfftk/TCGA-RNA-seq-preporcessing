import os
import argparse
import time
from utils import run_command

def main():
    parser = argparse.ArgumentParser(
        description="""
Perform STAR first-pass RNA-seq alignment.

This script performs the initial STAR alignment using raw or trimmed FASTQ files.
It supports single-end (SE) and paired-end (PE) reads, with multi-threading.

Workflow:
1. Reads accession numbers from input/SRRAccList.txt
2. Aligns reads to the specified genome using STAR first-pass parameters
3. Writes alignment stats to temporary/star/first_pass/stats/
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

    # Define input and output directories
    input_dir = "temporary/trimmomatic/output" if args.trimmed else "temporary/fastq_dump"
    first_pass_dir = "temporary/star/first_pass"
    genome_dir = "temporary/star/index/GRCh38_d1_vd1_v36_firstpass"
    stats_dir = os.path.join(first_pass_dir, "stats")

    os.makedirs(first_pass_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    # Read SRR accession list
    srr_list = [line.strip() for line in open("input/SRRAccList.txt")]

    print("=== STEP 1: First-pass STAR alignment ===")
    for srr in srr_list:
        start_time = time.time()

        # Prepare read files input for STAR
        if args.ending == "pe":
            read_files_in = f"{os.path.join(input_dir, srr+'_1.fastq')} {os.path.join(input_dir, srr+'_2.fastq')}"
        else:
            read_files_in = os.path.join(input_dir, srr + ".fastq")

        output_prefix = os.path.join(first_pass_dir, f"{srr}_")

        # Build STAR command
        command = (
            f"STAR "
            f"--genomeDir {genome_dir} "
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
            f"--readFilesCommand cat "
            f"--outFilterMatchNminOverLread 0.33 "
            f"--outFilterScoreMinOverLread 0.33 "
            f"--sjdbOverhang 100 "
            f"--outSAMstrandField intronMotif "
            f"--outSAMtype None "
            f"--outSAMmode None "
            f"2> {stats_dir}/{srr}_firstpass_stats.txt"
        )

        run_command(command, args.verbose)
        print(f"Completed first-pass for {srr} in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()
import os
import argparse
import time
from utils import run_command

def main():
    parser = argparse.ArgumentParser(
        description="""
Perform STAR first-pass RNA-seq alignment.

This script performs the initial STAR alignment using raw or trimmed FASTQ files.
It supports single-end (SE) and paired-end (PE) reads, with multi-threading.

Workflow:
1. Reads accession numbers from input/SRRAccList.txt
2. Aligns reads to the specified genome using STAR first-pass parameters
3. Writes alignment stats to temporary/star/first_pass/stats/
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

    # Define input and output directories
    input_dir = "temporary/trimmomatic/output" if args.trimmed else "temporary/fastq_dump"
    first_pass_dir = "temporary/star/first_pass"
    genome_dir = "temporary/star/index/GRCh38_d1_vd1_v36_firstpass"
    stats_dir = os.path.join(first_pass_dir, "stats")

    os.makedirs(first_pass_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    # Read SRR accession list
    srr_list = [line.strip() for line in open("input/SRRAccList.txt")]

    print("=== STEP 1: First-pass STAR alignment ===")
    for srr in srr_list:
        start_time = time.time()

        # Prepare read files input for STAR
        if args.ending == "pe":
            read_files_in = f"{os.path.join(input_dir, srr+'_1.fastq')} {os.path.join(input_dir, srr+'_2.fastq')}"
        else:
            read_files_in = os.path.join(input_dir, srr + ".fastq")

        output_prefix = os.path.join(first_pass_dir, f"{srr}_")

        # Build STAR command
        command = (
            f"STAR "
            f"--genomeDir {genome_dir} "
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
            f"--readFilesCommand cat "
            f"--outFilterMatchNminOverLread 0.33 "
            f"--outFilterScoreMinOverLread 0.33 "
            f"--sjdbOverhang 100 "
            f"--outSAMstrandField intronMotif "
            f"--outSAMtype None "
            f"--outSAMmode None "
            f"2> {stats_dir}/{srr}_firstpass_stats.txt"
        )

        run_command(command, args.verbose)
        print(f"Completed first-pass for {srr} in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()


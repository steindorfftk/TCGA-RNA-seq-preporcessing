import os
import argparse
import subprocess
import time
from pathlib import Path

def run_command(command, verbose):
    """Run a shell command and print it if verbose is enabled."""
    if verbose:
        print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(
        description="""
Trim bases from sequencing reads using Trimmomatic and run FastQC on the trimmed reads.

This script supports both single-end (SE) and paired-end (PE) reads. 
You can trim bases from either the start or the end of reads, optionally 
specify minimum read length, maximum cropping length, and the number of threads.

Workflow:
1. Reads accession numbers from input/SRRAccList.txt
2. Trims reads according to user-specified parameters
3. Runs FastQC on trimmed reads
"""
    )

    parser.add_argument(
        "-n", "--num_bases", type=int, required=True,
        help="Number of bases to remove."
    )
    parser.add_argument(
        "-p", "--position", choices=["start", "end"], required=True,
        help="Remove bases from the start or end of the reads."
    )
    parser.add_argument(
        "-e", "--ending", choices=["se", "pe"], required=True,
        help="Specify if the reads are single-end (SE) or paired-end (PE)."
    )
    parser.add_argument(
        "-nv", "--verbose", action="store_true", default=True,
        help="Enable verbose output (default: True)."
    )
    parser.add_argument(
        "-max", "--max_length", type=int,
        help="Maximum length to crop (used when trimming from end)."
    )
    parser.add_argument(
        "-th", "--threads", type=int, default=6,
        help="Number of threads to use (default: 6)."
    )
    parser.add_argument(
        "-min", "--minlen", type=int,
        help="Minimum length of reads to keep after trimming (optional)."
    )

    args = parser.parse_args()

    trimmomatic_path = "temporary/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
    
    # Determine Trimmomatic crop option
    trimmomatic_option = "HEADCROP" if args.position == "start" else "CROP"

    # Adjust num_bases if trimming from end with a maximum length
    if args.position == 'end' and args.max_length:
        args.num_bases = args.max_length - args.num_bases

    # Read SRR accession list
    srr_list = [line.strip() for line in open('input/SRRAccList.txt', 'r')]

    # Create output directories
    output_dir = "temporary/trimmomatic/output"
    os.makedirs(output_dir, exist_ok=True)
    fastqc_output = "temporary/fastqc/post_trimming_fastqc"
    os.makedirs(fastqc_output, exist_ok=True)

    for srr in srr_list:
        input_prefix = f'temporary/fastq_dump/{srr}'
        output_prefix = f'{output_dir}/{srr}'

        start_time = time.time()

        # Build Trimmomatic command
        if args.ending == "pe":
            trimmomatic_command = (
                f'java -jar {trimmomatic_path} PE '
                f'{input_prefix}_1.fastq {input_prefix}_2.fastq '
                f'{output_prefix}_1.fastq {output_prefix}_1_unpaired.fastq '
                f'{output_prefix}_2.fastq {output_prefix}_2_unpaired.fastq '
                f'{trimmomatic_option}:{args.num_bases} '
                f'-threads {args.threads}'
            )
        else:  # SE
            trimmomatic_command = (
                f'java -jar {trimmomatic_path} SE '
                f'{input_prefix}.fastq '
                f'{output_prefix}.fastq '
                f'{trimmomatic_option}:{args.num_bases} '
                f'-threads {args.threads}'
            )

        # Add MINLEN if provided
        if args.minlen:
            trimmomatic_command += f' MINLEN:{args.minlen}'

        # Run trimming
        run_command(trimmomatic_command, args.verbose)

        # Time taken for trimming
        elapsed_time = time.time() - start_time
        if args.verbose:
            print(f"Time taken for {srr}: {elapsed_time:.2f} seconds")

        # Run FastQC on trimmed files
        fastq_files = list(Path(output_dir).glob(f"{srr}*.fastq"))
        for fastq_file in fastq_files:
            if fastq_file.stat().st_size == 0:
                if args.verbose:
                    print(f"File {fastq_file.name} is empty. Deleting it.")
                fastq_file.unlink()
            else:
                if args.verbose:
                    print(f"Running FastQC on {fastq_file.name}...")
                run_command(f"fastqc -t {args.threads} -o {fastqc_output} {fastq_file}", args.verbose)

    if args.verbose:
        print("FastQC analysis completed.")

if __name__ == "__main__":
    main()


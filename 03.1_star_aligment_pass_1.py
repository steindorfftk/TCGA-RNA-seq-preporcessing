import os
import argparse
import time
from utils import run_command

def main():
    parser = argparse.ArgumentParser(description="STAR first-pass RNA-seq alignment.")
    parser.add_argument("-o", "--organism", choices=["Homo", "Mus"], required=True)
    parser.add_argument("-e", "--ending", choices=["se", "pe"], required=True)
    parser.add_argument("-nv", "--verbose", action="store_false", help="Disable verbose output.")
    parser.add_argument("-t", "--trimmed", action="store_true", help="Use trimmed FASTQs.")
    parser.add_argument("-th", "--threads", type=int, default=6)
    args = parser.parse_args()

    input_dir = "temporary/trimmomatic/output" if args.trimmed else "temporary/fastq_dump"
    first_pass_dir = "temporary/star/first_pass"
    genome_dir = "temporary/star/index/GRCh38_d1_vd1_v36_firstpass"
    stats_dir = os.path.join(first_pass_dir, "stats")

    os.makedirs(first_pass_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    srr_list = [line.strip() for line in open("input/SRRAccList.txt")]

    print("=== STEP 1: First-pass STAR alignment ===")
    for srr in srr_list:
        start_time = time.time()
        if args.ending == "pe":
            read_files_in = f"{os.path.join(input_dir, srr+'_1.fastq')} {os.path.join(input_dir, srr+'_2.fastq')}"
        else:
            read_files_in = os.path.join(input_dir, srr + ".fastq")

        output_prefix = os.path.join(first_pass_dir, f"{srr}_")
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


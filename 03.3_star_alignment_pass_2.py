import os
import argparse
import time
from utils import run_command

def main():
    parser = argparse.ArgumentParser(description="STAR second-pass RNA-seq alignment.")
    parser.add_argument("-o", "--organism", choices=["Homo", "Mus"], required=True)
    parser.add_argument("-e", "--ending", choices=["se", "pe"], required=True)
    parser.add_argument("-nv", "--verbose", action="store_false", help="Disable verbose output.")
    parser.add_argument("-t", "--trimmed", action="store_true", help="Use trimmed FASTQs.")
    parser.add_argument("-th", "--threads", type=int, default=6)
    args = parser.parse_args()

    input_dir = "temporary/trimmomatic/output" if args.trimmed else "temporary/fastq_dump"
    second_pass_index_dir = f"temporary/star/index/GRCh38_d1_vd1_v36_secondpass"
    second_pass_output_dir = "temporary/star/second_pass_output"
    stats_dir = os.path.join(second_pass_output_dir, "stats")
    
    os.makedirs(second_pass_output_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    srr_list = [line.strip() for line in open("input/SRRAccList.txt")]

    print("=== STEP 4: Second-pass STAR alignment ===")
    for srr in srr_list:
        start_time = time.time()
        if args.ending == "pe":
            read_files_in = f"{os.path.join(input_dir, srr+'_1.fastq')} {os.path.join(input_dir, srr+'_2.fastq')}"
        else:
            read_files_in = os.path.join(input_dir, srr + ".fastq")

        output_prefix = os.path.join(second_pass_output_dir, f"{srr}_")
        rg_line = f"ID:{srr} PL:ILLUMINA LB:{args.organism} SM:{srr}"

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

        run_command(command, args.verbose)
        print(f"Completed second-pass for {srr} in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()


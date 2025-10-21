import os
import argparse
from glob import glob
from utils import run_command

def main():
    parser = argparse.ArgumentParser(description="Generate STAR second-pass genome index.")
    parser.add_argument("-o", "--organism", choices=["Homo", "Mus"], required=True)
    parser.add_argument("-th", "--threads", type=int, default=6)
    parser.add_argument("-nv", "--verbose", action="store_false", help="Disable verbose output.")
    args = parser.parse_args()

    first_pass_dir = "temporary/star/first_pass"
    second_pass_index_dir = f"temporary/star/index/GRCh38_d1_vd1_v36_secondpass"
    reference_fasta = f"temporary/star/genome/GRCh38.d1.vd1.fa"

    os.makedirs(second_pass_index_dir, exist_ok=True)

    print("=== STEP 2: Collect SJ.out.tab files ===")
    sj_files = glob(os.path.join(first_pass_dir, "*SJ.out.tab"))
    if not sj_files:
        raise FileNotFoundError("No SJ.out.tab files found. Run first-pass first.")

    print(f"Found {len(sj_files)} SJ.out.tab files.")
    print("\n=== STEP 3: Generate second-pass genome index ===")

    command = (
        f"STAR "
        f"--runMode genomeGenerate "
        f"--genomeDir {second_pass_index_dir} "
        f"--genomeFastaFiles {reference_fasta} "
        f"--sjdbOverhang 100 "
        f"--runThreadN {args.threads} "
        f"--sjdbFileChrStartEnd {' '.join(sj_files)}"
    )

    run_command(command, args.verbose)
    print("Second-pass genome index generated successfully.")

if __name__ == "__main__":
    main()


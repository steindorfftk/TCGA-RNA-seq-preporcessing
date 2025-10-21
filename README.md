# RNA-seq STAR & FeatureCounts Pipeline

This pipeline performs RNA-seq analysis starting from SRA accession numbers. It includes downloading SRA data, converting to FASTQ, quality control, trimming, STAR alignment (first and second pass), and featureCounts quantification.

---

## Overview

The pipeline performs the following steps:

1. **Download SRA data** from a list of accession numbers.  
2. **Convert SRA to FASTQ** files using `fasterq-dump`.  
3. **Quality Control (FastQC)** on raw FASTQ files.  
4. **Trim reads** using Trimmomatic (optional).  
5. **STAR first-pass alignment** to generate SJ.out.tab files.  
6. **STAR second-pass genome index generation** using SJ.out.tab files.  
7. **STAR second-pass alignment** with improved accuracy.  
8. **Feature counting** using `featureCounts`.

---

## Requirements

- Python 3.8+  
- [SRA Toolkit](https://github.com/ncbi/sra-tools)  
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
- [STAR](https://github.com/alexdobin/STAR)  
- [Subread / featureCounts](http://subread.sourceforge.net/)  

## Installation

Clone the repository:

```bash
git clone https://github.com/steindorfftk/TCGA-RNA-seq_preprocessing.git
cd TCGA-RNA-seq_preprocessing
```

Place your **SRR accession list** in:

```
input/SRRAccList.txt
```

Example format:

```
SRR123456
SRR123457
SRR123458
```

---

## Directory Structure

The pipeline generates the following directories:

```
input/                     # Input SRR list
temporary/
├─ fastq_dump/             # Converted FASTQ files from SRA
├─ fastqc/                 # FastQC reports
├─ trimmomatic/             # Trimmed reads
├─ star/
│  ├─ first_pass/          # First-pass STAR alignment
│  ├─ second_pass_output/  # Second-pass STAR BAM files
│  └─ index/               # STAR genome indexes (first and second pass)
output/
├─ *_featurecounts.txt     # Raw featureCounts output
├─ *_gene_counts.tabular   # Cleaned gene count tables
```

---

## Pipeline Steps

### 1. Download SRA & Convert to FASTQ

```bash
python 01_download.py -o Mus -e pe -m -nv
```

Options:

- `-o`: Organism (`Mus` or `Homo`)  
- `-e`: Sequencing type (`se` or `pe`)  
- `-m`: Delete SRA files after conversion to save memory  
- `-nv`: Verbose output  

Output:

- `temporary/fastq_dump/` – FASTQ files  
- `temporary/fastqc/download_fastqc/` – FastQC reports  

---

### 2. Trim Reads (Optional)

```bash
python 02_trim_reads.py -n 10 -p start -e pe -th 8 -min 20 -max 100
```

Options:

- `-n`: Number of bases to trim  
- `-p`: Position to trim (`start` or `end`)  
- `-e`: Sequencing type (`se` or `pe`)  
- `-th`: Threads  
- `-min`: Minimum read length to keep  
- `-max`: Maximum length for trimming  

Output:

- `temporary/trimmomatic/output/` – Trimmed FASTQ files  
- `temporary/fastqc/post_trimming_fastqc/` – FastQC reports  

---

### 3. STAR First-Pass Alignment

```bash
python 03_star_firstpass.py -o Homo -e pe -t -th 6
```

Options:

- `-o`: Organism (`Homo` or `Mus`)  
- `-e`: Sequencing type (`se` or `pe`)  
- `-t`: Use trimmed FASTQs  
- `-th`: Threads  

Output:

- `temporary/star/first_pass/` – STAR output files  
- `temporary/star/first_pass/stats/` – Alignment stats  

---

### 4. STAR Second-Pass Genome Index

```bash
python 04_star_secondpass_index.py -o Homo -th 6
```

Options:

- `-o`: Organism  
- `-th`: Threads  

Output:

- `temporary/star/index/GRCh38_d1_vd1_v36_secondpass/` – STAR genome index  

---

### 5. STAR Second-Pass Alignment

```bash
python 05_star_secondpass_alignment.py -o Homo -e pe -t -th 6
```

Options:

- `-o`: Organism  
- `-e`: Sequencing type  
- `-t`: Use trimmed FASTQs  
- `-th`: Threads  

Output:

- `temporary/star/second_pass_output/` – Sorted BAM files  
- `temporary/star/second_pass_output/stats/` – Alignment stats  

---

### 6. Feature Counting

```bash
python 06_featurecounts.py -e pe -lm -nv
```

Options:

- `-e`: End type (`se` or `pe`)  
- `-lm`: Delete BAM files after processing  
- `-nv`: Verbose output  

Output:

- `output/` – Cleaned gene count tables (`*_gene_counts.tabular`)  

---

## Notes

- Ensure all tools are installed and available in PATH.  
- Input FASTQs/SRR lists must match the organism reference genome.  
- Low-memory mode (`-m` or `-lm`) is recommended for large datasets.  
- Thread count can be adjusted using `-th` for faster performance.  

---

## Author

Thiago Steindorff – GitHub: [steindorfftk](https://github.com/steindorfftk)  

---

Program Installation

    1 - Download the program files: Click on '<> Code' and then on 'Download ZIP'
    2 - Place the 'DEG-pre-processing-main.zip' file on the installation directory of your choice.
    3 - Change into the installation directory and unzip the program files (unzip DEG-pre-processing-main.zip)

Required packages Installation
SRA Tool kit

    1 - Download SRA tool kit from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit. For Linux, click on 'Ubuntu Linux 64 bit architecture' to download the file.
    2 - Place the 'sratoolkit.X.X.X-ubuntu64.tar.gz' file on the installation directory of your choice. (Replace X.X.X for the appropriate downloaded version)
    3 - Change into the installation directory and uzip the program files (tar -xzf sratoolkit.X.X.X-ubuntu64.tar.gz)
    4 - Change into sratoolkit.X.X.X-ubuntu64 directory.
    5 - Add SRA Tool kit binaries file to path (export PATH=$PATH:$PWD/bin).
    6 - Test that it was added to path: run 'which fastq-dump'. It should have a similar output as: '/Users/JoeUser/sratoolkit.3.0.0-mac64/bin/fastq-dump'
    7 - Delete 'init.py' from 'DEG-pre-processing-main/temporary/sratoolkit/'
    8.1 - Run the command 'vdb-config -i'. Use tab- and space/enter keys to navigate and select.
    8.2 - Enable Remote Access in the main menu
    8.3 - In the cache tab, enable local file-caching
    8.4 - In the cache tab, add the path to '/DEG-pre-processing-main/temporary/sratoolkit' to the "Location of user-repository"
    8.5 - Save and exit
    9.1 - Test that toolkit is functional: run 'fastq-dump --stdout -X 2 SRR390728'.
    9.2 - The output must be exaclty this: 'Read 2 spots for SRR390728 Written 2 spots for SRR390728 @SRR390728.1 1 length=72 CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC +SRR390728.1 1 length=72 ;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&( @SRR390728.2 2 length=72 AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT +SRR390728.2 2 length=72 ;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262'

STAR (with conda)

    1 - Run 'conda install -c bioconda star'
    2 - Download the reference genomes and index from https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/  
    3 - Place "gencode.v36.annotation.gtf" in /temporary/star/annotation, "GRCh38.d1.vd1.fa" in /temporary/star/genome and the index as "GRCh38_d1_vd1_v36_firstpass" in /temporary/star/index.

featureCounts (with conda)

    1.1 - Make a directory to install miniconda in the installation directory of your choice (mkdir -p miniconda3)
    1.2 - Change into miniconda3 directory
    1.3 - Download miniconda by executing (wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh)
    1.4 - Run miniconda installation script (bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3)
    1.5 - Remove miniconda installation script (rm -rf ~/miniconda3/miniconda.sh)
    1.6 - Initiallize miniconda3 ($PWD/bin/conda init bash)
    2.1 - Download features counts (conda install -c bioconda subread)
    2.2 - Test that featureCounts was installed (which featureCounts). It should have a similar output as: '/home/User/miniforge3/bin/featureCounts'
    The gene annotations at 'DEG-pre-processing-main/feature_counts/data' are updated as of february 2024. They were downloaded from Rsubread repository: https://code.bioconductor.org/browse/Rsubread/tree/RELEASE_3_9/inst/annot/(mkdir -p miniconda3)
    1.2 - Change into miniconda3 directory
    1.3 - Download miniconda (wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh)
    1.4 - Run miniconda installation script (bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3)
    1.5 - Remove miniconda installation script (rm -rf ~/miniconda3/miniconda.sh)
    1.6 - Initiallize miniconda3 ($PWD/bin/conda init bash)
    2.1 - Download features counts (conda install -c bioconda subread)
    2.2 - Test that featureCounts was installed correctly by checking its location using the command: (which featureCounts). It should have a similar output as: '/home/User/miniforge3/bin/featureCounts'
    The gene annotations located at 'DEG-pre-processing-main/feature_counts/data' have been updated as of February 2024. They were retrieved from the Rsubread repository: https://code.bioconductor.org/browse/Rsubread/tree/RELEASE_3_9/inst/annot/

fastqc

    1 - Run 'sudo apt-get install fastqc'

samtools

    1 - Run 'sudo apt-get install samtools'

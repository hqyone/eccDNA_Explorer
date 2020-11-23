#!/bin/bash
# Check the evironment and install softwares

wdir=$(pwd);
if [[ ! -d "$wdir" ]]; then $wdir="$PWD"; fi
binDir=$wdir/bin
if [[ ! -d "$binDir" ]]; then mkdir $binDir; fi

# https://cutadapt.readthedocs.io/en/stable/installation.html
if [[ -x "$(command -v cutadapt)" ]]
then
    Cutadapt=`which cutadapt`
else
    Cutadapt="$HOME/.local/bin/cutadapt"
    if [[ -x "$(command -v $Cutadapt)" ]] || [[ -x "$(command -v cutadapt)" ]]
    then
        echo "Cutadapt is found at $Cutadapt"
    else
        echo "install Cutadapt"
        python3 -m pip install --user --upgrade cutadapt
    fi
fi

# https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/https-cloning-errors
if [[ -x "$(command -v SeqPrep)" ]]
then
    SeqPrep=`which SeqPrep`
else
    SeqPrep="$binDir/SeqPrep/SeqPrep"
    if [[ -x "$(command -v $SeqPrep)" ]]
    then
        echo "SeqPrep is found at $SeqPrep"
    else
        echo "install SeqPrep ... "
        bin_dir="$binDir"
        if [[ -d $bin_dir ]]
        then
            echo "$bin_dir exist skip create the dir"
        else
            mkdir $bin_dir
        fi
        cd $bin_dir
        git clone https://github.com/jstjohn/SeqPrep.git
        cd "$bin_dir/SeqPrep"
        make
        echo "SeqPrep was installed at $binDir/SeqPrep/SeqPrep"
        cd $wdir
    fi
fi

# FastUniq : https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz
if [[ -x "$(command -v fastuniq)" ]]
then
    FastUniq=`which fastuniq`
else
    FastUniq="$binDir/FastUniq/source/fastuniq"
    if [[ -x "$(command -v $FastUniq)" ]]
    then
        echo "FastUniq is found at $FastUniq"
    else
        echo "install FastUniq ... "
        bin_dir="$binDir"
        if [[ -d $bin_dir ]]
        then
            echo "$bin_dir exist skip create the dir"
        else
            mkdir $bin_dir
        fi
        cd $bin_dir
        wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz
        tar -xzvf FastUniq-1.1.tar.gz
        rm FastUniq-1.1.tar.gz
        cd "$bin_dir/FastUniq/source"
        make
        echo "FastUniq was installed at $bin_dir/FastUniq/source/fastuniq"
        cd $wdir
    fi
fi


# BWA
if [[ -x "$(command -v bwa)" ]]
then
    BWA=`which bwa`
    echo "bwa is found $BWA"
else
    BWA="$binDir/bwa/bwa"
    if [[ -x "$(command -v $BWA)" ]] || [[ -x "$(command -v bwa)" ]]
    then
        echo "bwa is found $BWA"
    else
        echo "install bwa ... "
        bin_dir="$binDir"
        if [[ -d $bin_dir ]]
        then
            echo "$bin_dir exist skip create the dir"
        else
            mkdir $bin_dir
        fi
        cd $bin_dir
        git clone https://github.com/lh3/bwa.git
        cd "$bin_dir/bwa"
        make
        echo "bwa was installed at $bin_dir/bwa"
        cd $wdir
    fi
fi


# BedTools
if [[ -x "$(command -v bedtools)" ]]
then
    BedTools=`which bedtools`
    echo "BedTools is found at $BedTools"
else
    BedTools="$binDir/bedtools2/bin/bedtools"
    if [[ -x "$(command -v $BedTools)" ]] || [[ -x "$(command -v samtools)" ]]
    then
        echo "BedTools is found at $BedTools"
    else
        echo "install BedTools ... "
        bin_dir="$binDir"
        if [[ -d $bin_dir ]]
        then
            echo "$bin_dir exist skip create the dir"
        else
            mkdir $bin_dir
        fi
        cd $bin_dir
        git clone https://github.com/arq5x/bedtools2.git
        cd "$bin_dir/bedtools2"
        make
        echo "BedTools was installed at $bin_dir/bedtools2/bin"
        cd $wdir
    fi
fi

# Samtools
if [[ -x "$(command -v samtools)" ]]
then
    SAMTools=`which samtools`
    echo "SAMTools is found at $SAMTools"
else
    SAMTools="$binDir/samtools/samtools"
    if [[ -x "$(command -v $SAMTools)" ]]
    then
        echo "SAMTools is found at $SAMTools"
    else
        echo "install SAMTools ... "
        bin_dir="$binDir"
        if [[ -d $bin_dir ]]
        then
            echo "$bin_dir exist skip create the dir"
        else
            mkdir $bin_dir
        fi
        cd $bin_dir
        git clone https://github.com/samtools/htslib.git
        cd "$bin_dir/htslib"
        make
        make install
        cd $bin_dir
        git clone https://github.com/samtools/samtools.git
        cd "$bin_dir/samtools"
        make
        make install
        echo "SAMTools was installed at $bin_dir/samtools/samtools"
        cd $wdir
    fi
fi

# Human Genome FASTA and GTF downloading and indexing
cd $wdir
#genome_dir="$wdir/genome"
genome_dir="/home/hqyone/mnt/3t/rna_seq/genome"
gfa_url="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
gfa_gz_file=`basename $gfa_url`
GenomeFASTA="$genome_dir/${gfa_gz_file/.fa.gz/.fa}"

gtf_url="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz"
gtf_gz_file=`basename $gtf_url`
GenomeGTF="$genome_dir/${gtf_gz_file/.gtf.gz/.gtf}"

if ! [[ -d $genome_dir ]];then
    mkdir "$genome_dir"
fi

cd $genome_dir
if [[ -f $GenomeFASTA ]];then
    echo "Find the genome fasta file at $GenomeFASTA"
else
    wget $gfa_url -O $gfa_gz_file
    gunzip -c $gfa_gz_file > $GenomeFASTA
    rm $gfa_gz_file
fi

if [[ -f $GenomeGTF ]];then
    echo "Find the genome GTF file at $GenomeGTF"
else
    wget $gtf_url -O $gtf_gz_file
    gunzip -c $gtf_gz_file > $GenomeGTF
    rm $gtf_gz_file
fi

# Index genome by BWA
#bwa_index_file=${GenomeFASTA%".fa"}.bwt
bwa_index_file="$GenomeFASTA.bwt"
echo $bwa_index_file
if ! [ -f $bwa_index_file ];then
    echo "Indexing genome by bwa"
    eval "$BWA index $GenomeFASTA"
    #echo $output
fi

# Index genome by SAMTOOLS
#samtools_index_file=${GenomeFASTA%".fa"}.bai
samtools_index_file="$GenomeFASTA.fai"
echo $samtools_index_file
if ! [ -f $samtools_index_file ];then
    echo "Indexing genome by SAMTools"
    eval "$SAMTools faidx $GenomeFASTA"
    #echo $output
fi



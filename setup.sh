#!/bin/bash
# Check the evironment and install softwares

wdir=$(pwd);
if [[ ! -d "$wdir" ]]; then $wdir="$PWD"; fi
binDir=$wdir/bin
if [[ ! -d "$binDir" ]]; then mkdir $binDir; fi

# https://cutadapt.readthedocs.io/en/stable/installation.html
Cutadapt="$HOME/.local/bin/cutadapt"
if [[ -x "$(command -v $Cutadapt)" ]] || [[ -x "$(command -v cutadapt)" ]]
then
    echo "Cutadapt is found at $Cutadapt"
else
    echo "install Cutadapt"
    python3 -m pip install --user --upgrade cutadapt
fi

# https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/https-cloning-errors
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

# BWA
BWA="$binDir/bwa/bwa"
if [[ -x "$(command -v $BWA)" ]]
then
    echo "bwa is found at $BWA"
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


# BedTools
BedTools="$binDir/"
if [[ -x "$(command -v $BedTools)" ]]
then
    echo "BedTools is found at $BedTools"
else
    echo "install BedTools ... "
    bin_dir="$binDir/bedtools2"
    if [[ -d $bin_dir ]]
    then
        echo "$bin_dir exist skip create the dir"
    else
        mkdir $bin_dir
    fi
    cd $bin_dir
    wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
    tar -zxvf bedtools-2.29.2.tar.gz
    cd "$bin_dir/bedtools2"
    make
    echo "BedTools was installed at $bin_dir/bedtools2"
    cd $wdir
    
    echo "sss"
    echo "ooo"
fi
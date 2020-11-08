#!/bin/bash

wdir=$(pwd);
if [[ ! -d "$wdir" ]]; then $wdir="$PWD"; fi
. "$wdir/setup.sh"

echo " #### Cutadapt : $Cutadapt"
echo " #### SeqPrep : $SeqPrep"
echo " #### BWA : $BWA"

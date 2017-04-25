#!/bin/bash
#$ -l mem=4G,time=4:: -N paraDemo -S /bin/bash -cwd 

# This script is a demo shell script with parameters
# add file readme information here
# 
#    FqDir - (required) - A directory containing data files
#    Type - (required) - "csv" for comma seperate, "tsv" for tab seperate
#    OutNam - (optional) - A name for the output file. 
#    Help - H - (flag) - get usage information


#list of required tools:
# standard linux library
###############################################################

usage="
<-t 1-2>* para_sh.sh -i <InputDirectory> -t <Type> -o <OutputName>
     -i (required) - Path to directory containing fastq.gz files
     -t (required) - \"csv\" for comma \"tsv\" for tab seperate files
     -o (optional) - Output filename - if not provided the directory name will be used
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:t:o:H opt; do
    case "$opt" in
        i) FqDir="$OPTARG";;
        t) Type="$OPTARG";; 
        o) OutNam="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#check directory exists
if [[ ! -d $FqDir ]]; then
    echo  "Need provide a directory"
    echo $usage
    exit
fi

#check for csv/tsv specification
if [[ "$Type" != "PE" ]] && [[ "$Type" != "SE" ]]; then
    echo  "Need to specify paired-end or single-end"
    echo $usage
    exit
fi


### do something you want

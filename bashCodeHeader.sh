
#!/bin/bash
# This script generates collapsed IGV plots given variant position and bam location 
# the plots will be in the same folder of the given variant file
# example: sh Generate_IGV_plots.sh -b example/bam.txt -v example/IGV_variants.txt
#

###############################################################
#set default arguments
usage="
     -v (required) - Table containing the path to the fastq file and the RG read header
     -b (required) - shell file containing variables with locations of reference files and resource directories (WES_Pipeline_References.b37.sh)
     -E (flag) - generate IGV plots in expand mode
     -H (flag) - echo this message and exit
"

#get arguments
ExpandMode="false"
while getopts v:b:EH opt; do
    case "$opt" in
        v) INDELS="$OPTARG";;
        b) BAM="$OPTARG";;
                E) ExpandMode="true";;
        H) echo "$usage"; exit;;
    esac
done

echo $ExpandMode


### input information
INDELS=`readlink -f $INDELS`
BAM=`readlink -f $BAM`
BamNam=$(echo $INDELS | sed s/.txt// )
DIR=$BamNam".IGV_plots/" ## figure output folder
mkdir -p $DIR


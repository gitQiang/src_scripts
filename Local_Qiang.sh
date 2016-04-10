#!/bin/bash

BAMS="/Volumes/TwoT/Desktop/breastcancer/TCGA/bambatch3.txt"
#BAMS="/Volumes/TwoT/Desktop/breastcancer/TCGA/testbam.txt"
RefFil="/Volumes/TwoT/Desktop/breastcancer/TCGA/Exome-vc-pipelineV2/WES_Pipeline_References.b37plusDecoy_Qiang_local.sh"
TgtBed="/Volumes/TwoT/Desktop/breastcancer/TCGA/CBL_target.bed"

while read line
do
	## step 1 run ExmAln.2.HaplotypeCaller_GVCFmode.sh to generate gvcf files
	./ExmAln.2.HaplotypeCaller_GVCFmode.sh -i $line -r $RefFil -t $TgtBed
	
	## step 2 run ExmVC.1hc.GenotypeGVCFs.sh to get VCF files
	subject=`basename $line | sed 's/.bam$//'`
	#echo $subject
	./ExmVC.1hc.GenotypeGVCFs.sh -i $subject.g.vcf.gz -r $RefFil -t $TgtBed -X
	
	## Step 3 run ExmVC.3.AnnotateVCF.sh to annotate the variants
	cp $subject.g.vcf.gz.splitfiles/$subject.g.vcf.gz.1.vcf ./
	#./ExmVC.3.AnnotateVCF.sh -i $subject.g.vcf.gz.1.vcf -r $RefFil -X 
	
	## Step 4 run ExmVC.4.RecalibrateVariantQuality.sh to variant Quality
	#./ExmVC.4.RecalibrateVariantQuality.sh -i $subject.g.vcf.gz.1.vcf -r $RefFil -X
	
done < "$BAMS"










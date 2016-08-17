#! /bin/bash

#$ -N sample

#$ -q all.q

###$ -wd

#$ -cwd

#$ -e $JOB_ID.$JOB_NAME.e

#$ -o $JOB_ID.$JOB_NAME.o

#$ -pe smp 4

#$ -V

name=sample



#bwa aln_1

bwa aln -t $NSLOTS -f $name\_1.sai /shared/data/GATKbundle/hg19.fasta $name\_1.fastq.gz



#bwa aln_2

bwa aln -t $NSLOTS -f $name\_2.sai /shared/data/GATKbundle/hg19.fasta $name\_2.fastq.gz



#bwa sampe

bwa sampe -f $name\.sam /shared/data/GATKbundle/hg19.fasta $name\_1.sai $name\_2.sai $name\_1.fastq.gz $name\_2.fastq.gz



##picard SortSam

picard-tools SortSam

SO=coordinate I=$name\.sam O=$name\.bam

VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true



#picard AddOrReplaceReadGroups

picard-tools AddOrReplaceReadGroups

I=$name\.bam O=$name\_RG.bam

SORT_ORDER=coordinate

VALIDATION_STRINGENCY=SILENT

RGID=$name RGLB=$name RGPL=Illumina RGPU=$name RGSM=$name



#picard MarkDuplicates

picard-tools MarkDuplicates

I=$name\_RG.bam O=$name\_dupRemoved.bam

METRICS_FILE=$name\_dupRemoved.met

VALIDATION_STRINGENCY=SILENT

REMOVE_DUPLICATES=true



#samtools

samtools index $name\_dupRemoved.bam



#GATK RealignerTargetCreator

java -Xmx4g -jar /shared/tools/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar 

	-R /shared/data/GATKbundle/hg19.fasta 

	-T RealignerTargetCreator

	-known /shared/data/GATKbundle/1000G_phase1.indels.hg19.new.vcf

	-known /shared/data/GATKbundle/Mills_and_1000G_gold_standard.indels.hg19.new.vcf

	-I $name\_dupRemoved.bam -o $name\_dupRemoved.intervals





#GATK IndelRealigner

java -Xmx4g -jar /shared/tools/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

	-R /shared/data/GATKbundle/hg19.fasta

	-T IndelRealigner

	-targetIntervals $name\_dupRemoved.intervals

	-I $name\_dupRemoved.bam

	-o $name\_dupRemoved_realigned.bam



#picard FixMateInformation

picard-tools FixMateInformation

INPUT=$name\_dupRemoved_realigned.bam

OUTPUT=$name\_dupRemoved_realigned_fixed.bam

SO=coordinate

VALIDATION_STRINGENCY=SILENT

CREATE_INDEX=true



#GATK BaseRecalibrator

java -Xmx4g -jar /shared/tools/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

	-R /shared/data/GATKbundle/hg19.fasta

	-T BaseRecalibrator

	-I $name\_dupRemoved_realigned_fixed.bam

	-knownSites /shared/data/GATKbundle/dbsnp_137.hg19.excluding_sites_after_129.new.vcf

	-o $name\_recal_data.grp



#GATK PrintReads 

java -Xmx4g -jar /shared/tools/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

	-R /shared/data/GATKbundle/hg19.fasta

	-T PrintReads

	-I $name\_dupRemoved_realigned_fixed.bam

	-BQSR $name\_recal_data.grp

	-o $name\_dupRemoved_realigned_fixed_recalib.bam


#$ -S /bin/bash

#$ -N sample

#$ -q all.q

###$ -wd

#$ -cwd

#$ -e $JOB_ID.$JOB_NAME.e

#$ -o $JOB_ID.$JOB_NAME.o

#$ -pe smp 4

#$ -V

#


#MuTect for SNP

java -Xmx4g -jar /shared/tools/muTect-1.1.7/mutect-1.1.7.jar

	--analysis_type MuTect

	--reference_sequence /shared/data/GATKbundle/ucsc.hg19.fasta

	--cosmic /shared/tools/muTect-1.1.7/b37_cosmic_v54_120711_ucsc.hg19.compatible.vcf

	--dbsnp /shared/tools/muTect-1.1.7/dbsnp_132_b37.leftAligned_ucsc.hg19.compatible.vcf

	--input_file:normal $name\N_dupRemoved_realigned_fixed_recalib.bam

	--input_file:tumor $name\T_dupRemoved_realigned_fixed_recalib.bam

	--out call.stats.$name\NT.out

	--coverage_file coverage.$name\NT.wig.txt

	--vcf $name\NT.snp.vcf



#Filtering

grep "PASS" $name\NT.snp.vcf > anno.$name\NT.snp.vcf



#ANNOVAR

/shared/tools/annovar/table_annovar.pl anno.$name\NT.snp.vcf

/shared/tools/annovar/humandb/

-buildver hg19

-out myanno.$name\NT.snp

-remove

-protocol refGene,cytoBand,genomicSuperDups,gwasCatalog,1000g2014oct_all,1000g2014oct_eas,

cosmic70,esp6500siv2_all,esp6500siv2_ea,exac03,snp138,cg46,clinvar_20150330,nci60,ljb26_all

-operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f

-nastring . 

-vcfinput 

-otherinfo



#SomaticIndelDetector for InDel

java -Xmx4g -jar /shared/data/GenomeAnalysisTK-2.1-13/GenomeAnalysisTK.jar

	-R /shared/data/GATKbundle/ucsc.hg19.fasta

	-T SomaticIndelDetector

	-o $name\NT.indel.vcf

	-verbose $name\NT.indels.txt

	-I:normal $name\N_dupRemoved_realigned_fixed_recalib.bam

	-I:tumor $name\T_dupRemoved_realigned_fixed_recalib.bam



#Filtering

grep "SOMATIC" $name\NT.indel.vcf > anno.$name\NT.indel.vcf



#ANNOVAR

shared/tools/annovar/table_annovar.pl anno.$name\NT.indel.vcf /shared/tools/annovar/humandb/

-buildver hg19

-out myanno.$name\NT.indel

-remove

-protocol refGene,cytoBand,genomicSuperDups,gwasCatalog,1000g2014oct_all,1000g2014oct_eas,

cosmic70,esp6500siv2_all,esp6500siv2_ea,exac03,snp138,cg46,clinvar_20150330,nci60,ljb26_all

-operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f

-nastring . 

-vcfinput 

-otherinfo

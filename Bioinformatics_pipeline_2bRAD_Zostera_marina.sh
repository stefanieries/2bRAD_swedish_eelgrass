# 2bRAD pipline log for Zostera marina
# Edited by Stefanie Ries Septmeber 2023
# Based on the last edited version by Ellika Faust November 2022

###############################################
### 	     	Rackham		            ###
###############################################
# user guides
# https://uppmax.uu.se/support-sv/user-guides/slurm-user-guide/
# https://uppmax.uu.se/support-sv/user-guides/rackham-user-guide/

#downloading data from dds client
dds data get -p snpseq00155  -a --verify-checksum -d VJ-3470

uquota #file system usage on Rackham
projinfo #CPU hrs
jobinfo -u sries
finishedjobinfo -j
squeue -u sries
projsummary[project id] #overview

#log in to Rackham
ssh -X sries@rackham.uppmax.uu.se

###########
## There's a problem using wildcards (*) for scp on a Mac. In short OsX thinks * means something else.
## You can easily fix this by putting quotation marks around the filepath.
## scp 'usr@rackham.uppmax.uu.se:/dir/path/*.gz' .
## More details here: https://unix.stackexchange.com/questions/27419/how-to-use-wildcards-when-copying-with-scp (bearbeitet) 
###########

cd /proj/snic2022-23-656/2bRAD_analysis/fastqs # my data

# info given by UPPMAX verify the integrity of the downloaded files with checksums
# for each checksum file, go to the folder where it is located and verify the checksums within
#mv checksums.md5 one folder up
cp checksums.md5 ..
md5sum -c checksums.md5

# -n nodes
# -t h:mm:ss
# -A project
# -p partitioning

interactive -n 1 -t 5:00:00 -A snic2022-22-1045 # Steffis

# for a job submission

#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p node -n 64
#SBATCH -t 2:00:00
#SBATCH -J jobname
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se

# to submit
# sbatch <filename>


#####################################
###     DOWNLOAD 2bRAD SCRIPTS    ###
#####################################

# downloading and installing all 2bRAD scripts in $HOME/bin (or change to whatever directory you want)
cd
mkdir bin
cd ~/bin
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
# move scripts to ~/bin from sub-directories
mv 2bRAD_denovo/* .
# remove now-empty directory
rm -rf 2bRAD_denovo
rm -rf 2bRAD_GATK

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# adding ~/bin to your $PATH
cd
nano .bashrc
# paste this where appropriate (note: .bashrc configuration might be specific to your cluster, consult your sysadmin if in doubt)
export PATH=$HOME/bin:$PATH

# Ctl-o, Ctl-x  (to save and exit in nano)
# log out and re-login to make sure .bashrc changes took effect

# does it work?
# try running a script from $HOME:
cd
2bRAD_trim_launch.pl
# if you get "command not found" something is wrong

#####################################
### DOWNLOAD READS FROM BASESPACE ###
#####################################

#organize fastqs
mkdir fastqs
cd fastqs

# copy them like this
rsync --progress /proj/snic2022-23-656/VJ-3470/files/VJ-3470/221209_A00181_0603_AHGYLJDMXY/221209_A00181_0603_AHGYLJDMXY/Sample_VJ-3470*/* /proj/snic2022-23-656/2bRAD_analysis/fastqs

ls ./*L001_R1_001.fastq.gz | wc -l

############################
####### INDEX GENOME #######
############################

# ran this interactively
module load bioinfo-tools
module load bowtie2
module load samtools
module load picard/2.23.4

export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.dict"
bowtie2-build $REFERENCE_GENOME $REFERENCE_GENOME
samtools faidx $REFERENCE_GENOME
java -jar $PICARD_ROOT/picard.jar CreateSequenceDictionary R=$REFERENCE_GENOME  O=$REFERENCE_GENOME_DICT


######################################
#### CONCATENTATE LANE DUPLICATES ####
######################################

# This is a special case for VF-3360 data where there was an issue with demultiplexing and we have 6 files per samples rather than 2.
#for file in ./*_barcode_Unknown_L00*_R1_001.fastq.gz
#do echo `cat ${file} ${file/_barcode_Unknown_/_CGGGCT_} ${file/_barcode_Unknown_/_AGATCT_} > ${file/_barcode_Unknown_/_}`
#done

>concat.sh
'#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core -n 1
#SBATCH -J concat
#SBATCH -o concat.%J.out
#SBATCH -e concat.%J.err
#SBATCH --time=01:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se

#forward
for file in ./*L001_R1_001.fastq.gz
do echo `cat ${file} ${file/_L001_/_L002_} > ${file/_L001_R1_001.fastq.gz/}_R1.fq.gz`
done'>>concat.sh
# THis is a special case for VF-3360 data where there was an issue with demultiplexing and we have 6 files per samples rather than 2.
#rm *_barcode_Unknown_*fq.gz
#rm *_CGGGCT_*fq.gz
#rm *_AGATCT_*fq.gz


#CHECK CONCATENATENATION RESULTS MAKE SENSE
ls *R1.fq.gz | wc -l

################################
### parallel gzip and gunzip ###
###############################
>unpigz.sh
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core -n 8
#SBATCH -J unpigz
#SBATCH -o unpigz.%J.out
#SBATCH -e unpigz.%J.err
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se

unpigz -p 8 -v -f /proj/snic2022-23-656/2bRAD_analysis/fastqs/*gz'>>unpigz.sh

>pigz.sh
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core -n 8
#SBATCH -J pigz
#SBATCH -o pigz.%J.out
#SBATCH -e pigz.%J.err
#SBATCH --time=02:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se

pigz -p 8 -v -f /proj/snic2022-23-656/2bRAD_analysis/fastqs/*fq
# wait
# pigz -p 8 -v -f /proj/snic2022-23-656/2bRAD_analysis/fastqs/*tr0
# wait
#pigz -p 8 -v -f /proj/snic2022-23-656/2bRAD_analysis/fastqs/*trim
# wait
#pigz -p 8 -v -f /proj/snic2022-23-656/2bRAD_analysis/fastqs/*sam
# wait
#pigz -p 8 -v -f /proj/snic2022-23-656/2bRAD_analysis/fastqs/*bam
'>>pigz.sh

######################################
###            RUN FASTQ           ###
######################################

#RUN FASTQC
module load bioinfo-tools
module load FastQC
mkdir Fastqc_Results/
> runFQC
echo '#!/bin/bash -l' > runFQC
for file in /proj/snic2022-23-656/2bRAD_analysis/fastqs/*.fq
do echo "fastqc -o /proj/snic2022-23-656/2bRAD_analysis/fastqs/Fastqc_Results -f fastq $file " >> runFQC
done

sh runFQC

ls Fastqc_Results/*zip | wc -l
ls Fastqc_Results/*html | wc -l

# copy fastQC results on local computer
scp sries@rackham.uppmax.uu.se:'/proj/pipelines_2023_data/nobackup/sries/QC_data/*zip' /Users/xriest/Documents/0_PhD_University_Gothenburg/2_PhD_Credits_Courses/2.1_Courses/2023_NMAR302_PopulationGeneticsCourse/bioinformatic_exercises/  

############################
###     PREP READS       ###
############################

# Trimming low-quality bases and remove pcr duplicates

# This method is sequential, far from optimal
# takes roughly 1 min and 30 sek per command
# ==========================
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 13:00:00
#SBATCH -J trimming
#SBATCH -o trimming.%J.out
#SBATCH -e trimming.%J.err
#SBTACH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se

module load bioinfo-tools
module load perl' > trimming.sh

trim2bRAD_dedup.pl fq >> trimming.sh

# do we have expected number of *.tr0 files created?
ll *R1.fq.tr0 | wc -l

# ===============================
> trimse.sh
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J trimse
#SBATCH -o trimse.%J.out
#SBATCH -e trimse.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se
module load bioinfo-tools
module load cutadapt' > trimse.sh

# for reference-based analysis: trimming poor quality bases off ends:
# removing reads with qualities at ends less than Q15
# mishas script says you can use -m 25 (minimum read length) if you are aligning to genome

for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse.sh;
done

# execute all commands in trimse file (serial or parallel using Launcher, if your system allows)
# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

#==============
# Mapping reads to reference (reads-derived fake one, or real) and formatting bam files

# for reference-based:
export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.dict"

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
> maps.sh
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J maps
#SBATCH -o maps.%J.out
#SBATCH -e maps.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se
module load bioinfo-tools
module load bowtie2' > maps.sh


2bRAD_bowtie2_launch.pl '\.trim$' $REFERENCE_GENOME >> maps.sh

# takes a veeeery long time, so I split them in seperat sbatch scripts, would really help to parallelize this
# remember to change the job names

head -n12 maps.sh > mapsVF.sh
grep 'VF' maps.sh >> mapsVF.sh
head -n12 maps.sh > mapsUJ.sh
grep 'UJ' maps.sh >> mapsUJ.sh

# execute all commands written to maps

cat mapshead.30843205.err mapstail.30843204.err > maps
ls *.sam | wc -l  # number should match number of trim files
grep 'unpaired' maps | wc -l # just make sure this is also the same as .sam files are written part by part

>alignmentRates
for F in `ls *trim`; do
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
# adding read groups, validating bams
#-------------------------------

# Add read groups to sam files before converting to bam
echo '#!/bin/bash -l '> addrg.sh # 1 min for 10 files
for i in *.trim.bt2.sam; do
	newfile="$(basename $i .fq.trim.bt2.sam)"
	echo "addReadGroup.pl -i $i -o ${newfile}.fq.trim.bt2.rg.sam -r ${newfile} -s ${newfile} -l ${newfile} -u ${i:0:7} -p Illumina -c NGI" >>addrg.sh;
done

sbatch -A snic2022-22-1045 -p core -n 1 -t 1:00:00 -o addrg.%J.out -e addrg.%J.err --mail-type ALL --mail-user stefanie.ries@gu.se addrg.sh

ls *bt2.sam | wc -l
ls *rg.sam | wc -l

# validate sam files
# this takes a long time so skip this step unless you run into an error

#>validateSams.sh # 3 min fr 10 samples
#for i in *rg.sam; do
#	echo "java -jar $PICARD_ROOT/picard.jar ValidateSamFile -I $i -MODE SUMMARY 1>>sam_validationsummary.out 2>>sam_validationsummary.err" >>validateSams.sh;
#done
#>sam_validationsummary.out
#>sam_validationsummary.err
# sh validateSams.sh


# convert, sort and index from sam to bam file
export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.dict"

echo '#!/bin/bash -l
module load bioinfo-tools
module load samtools'> s2b.sh # 2min 20 sec for 10 samples

for file in *rg.sam; do
echo "samtools sort --reference $REFERENCE_GENOME -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b.sh;
done

sbatch -A snic2022-22-971 -p core -n 1 -t 1:00:00 -o s2b.%J.out -e s2b.%J.err --mail-type ALL --mail-user stefanie.ries@gu.se s2b.sh

#################################

ls *rg.bam | wc -l

echo '#!/bin/bash -l
module load bioinfo-tools
module load picard'>validateBams.sh # 1 min 30 sec for 10 samples
>bam_validationsummary.out
>bam_validationsummary.err
for i in *rg.bam; do
	echo "java -jar \$PICARD_ROOT/picard.jar ValidateSamFile -I $i -R $REFERENCE_GENOME -MODE SUMMARY 1>>bam_validationsummary.out 2>>bam_validationsummary.err" >>validateBams.sh;
done

sbatch -A snic2022-22-971 -p core -n 1 -t 1:00:00 -o validateBams.%J.out -e validateBams.%J.err --mail-type ALL --mail-user stefaine.ries@gu.se validateBams.sh

cat bam_validationsummary.out

# If everything looks good, remove all sam files
# rm *bt2.sam
# rm *rg.sam

###########################################################################################################################################
# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them. #
###########################################################################################################################################

# Next we check the average depth to see if we proceed with GATK (>10x) or ANGSD (<10x)
# if your coverage is >10x, go to GATK section below

module load bioinfo-tools
module load samtools

> average_depth
for file in *.bam; do
echo "working on $file"
D=`samtools depth $file | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
echo "$file $D" >> average_depth;
done

# listing all bam filenames
ls *bam >bams
#----------- assessing base qualities and coverage depth

module load bioinfo-tools
module load samtools
module load ANGSD

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option)
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples (1053 samples x 10 = 10530)

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 10530 -minInd 10"

# T O   D O :
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

angsd -b bams -r Chr01 -GL 1 $FILTERS $TODO -P 1 -out dd

# scp dd to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs
# this took too long for over 1000 samples -> sent as job!

> basemapping_coveragedepth.sh
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00 # takes roughly 2,5 h fpr 1000 samples
#SBATCH -J basemapping_coveragedepth
#SBATCH -o basemapping_coveragedepth.%J.out
#SBATCH -e basemapping_coveragedepth.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se

module load bioinfo-tools
module load samtools
module load ANGSD

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 10530 -minInd 10"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

angsd -b bams -r Chr01 -GL 1 $FILTERS $TODO -P 1 -out dd

'>>basemapping_coveragedepth.sh

sbatch basemapping_coveragedepth.sh

# summarizing results (using modified script by Matteo Fumagalli)
module load RStudio
module load R_packages/4.1.1

Rscript ~/bin/plotQC.R prefix=dd

# copy dd files it in my mac folder

# proportion of sites covered at >5x:
cat quality.txt

#############################
#        G  A  T  K         #
#############################

#        G  A  T  K
# ("hard-call" genotyping, use only for high-coverage data, >10x after deduplication)

#-------------------------------------------------
# renaming individuals
#-------------------------------------------------
# NOTE: here I renamed some individual bam files that were mislabeled
# the folloing files were misslabelled and are here corrected with rename eg.
# rename YST-19 STE-18-rep *
# had to also add Readgroup again
# and index
# samtools index *bam

# Wrong name	Correct name
# YST-19 	    STE-18-rep
# BJO-4-Rep 	ALA-05-rep
# FUR-20-Rep 	HOR-11-rep
# GRO-8-Rep 	HOG-12-rep
# HOG-12-Rep 	GRO-08-rep
# HOR-11-Rep 	FUR-20-rep
# ALA-5-Rep 	BJO-04-rep
# STE-18-Rep 	YST-19
# -----------------------------------------------


export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.dict"

ls *.bam > bams

# writing command script with SLURM header (some fields might be different on your cluster, contact your IT people!)
>unigt.sh
echo '#!/bin/bash
#SBATCH -A snic2022-22-1045
#SBATCH -n 20
#SBATCH -o unigt.%J.out
#SBATCH -e unigt.%J.err
#SBATCH -t 48:00:00 # takes a long time!
#SBATCH --mail-type ALL
#SBATCH --mail-user stefanie.ries@gu.se
module load bioinfo-tools
module load GATK/3.8-0
java -Djava.io.tmpdir=$TMPDIR -jar $GATK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R /proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta -nct 20 \
--genotype_likelihoods_model SNP \' >unigt.sh
cat bams | perl -pe 's/(\S+\.bam)/-I $1 \\/' >> unigt.sh
echo '-o primary.vcf ' >> unigt.sh

sbatch unigt.sh

### in the following lines I added the number from my seascape dataset (all data -> incl. monitoring samples)
# Implementing a minDP3 filter here prior to VQSR
module load bioinfo-tools
module load vcftools
vcftools --vcf primary.vcf --minDP 3 --recode --recode-INFO-all --out primary_DP3
# After filtering, kept 1053 out of 1053 Individuals
# Outputting VCF file...
# After filtering, kept 63445 out of a possible 63445 Sites

#----------
# Variant quality score recalibration (VQSR)

# making a tab-delimited table of clone (replicate) sample pairs
#cat clonepairs.tab
cat clonepairs_oldnames.tab # called this "oldnames" to prevent confusion, when we uploaded the new file with renamed sampels
# paste names of bams that are clone pairs, tab delimited, one pair per line; for example
# Ctl-O , enter, Ctl-X

#----------
# special step: we needed to rename some clonpairs
# upload the clonepairs file with manually renamed names
scp /Users/xriest/Documents/clonepairs_renamed.tab sries@rackham.uppmax.uu.se:/proj/snic2021-23-661/2bRAD_analysis/fastqs/renamed_bams/
#----------

# extracting "true snps" subset (reproducible across replicates)
# parameter hetPairs can vary depending on replication scheme (3 is good when you have triplicates)

# had to run with default hetPairs =1 rather than 2 as we got to few variants and GATK threw an error
# also tried som other adjustments but none did the trick...
replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs_oldnames.tab altPairs=0 hetPairs=0  > vqsr.vcf
# 63445 total SNPs
# 8218 pass hets and match filters
# 2314 show non-reference alleles
# 8218 have alterantive alleles in at least 0 replicate pair(s)
# 8218 have matching heterozygotes in at least 0 replicate pair(s)
# 2296 polymorphic
# 2314 written

replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs_oldnames.tab hetPairs=1  > vqsr1.vcf
# 63445 total SNPs
# 8218 pass hets and match filters
# 2314 show non-reference alleles
# 2314 have alterantive alleles in at least 1 replicate pair(s)
# 2270 have matching heterozygotes in at least 1 replicate pair(s)
# 2270 polymorphic
# 2270 written

replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs_oldnames.tab hetPairs=2  > vqsr2.vcf
# 63445 total SNPs
# 8218 pass hets and match filters
# 2314 show non-reference alleles
# 2314 have alterantive alleles in at least 1 replicate pair(s)
# 621 have matching heterozygotes in at least 2 replicate pair(s)
# 621 polymorphic
# 621 written

# determining transition-transversion ratio for true snps (will need it for tranche calibration)
module load bioinfo-tools
module load vcftools
vcftools --vcf vqsr.vcf --TsTv-summary
# Ts/Tv ratio: 2.271
vcftools --vcf vqsr1.vcf --TsTv-summary
# Ts/Tv ratio: 2.801
vcftools --vcf vqsr2.vcf --TsTv-summary
# Ts/Tv ratio: 2.935

# creating recalibration models
export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
module load bioinfo-tools
module load GATK/3.8-0

# put your actual number into the next code chunk, --target_titv

# VariantRecalibrator = Build a recalibration model to score variant quality for filtering purposes
# This tool performs the first pass in a two-stage process called Variant Quality Score Recalibration (VQSR)

### OPTION: QD, DP, FS, MQ, SOR, MQRankSum
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar  -T VariantRecalibrator \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
-resource:repmatch,known=true,training=true,truth=true,prior=15  vqsr1.vcf \
-an QD -an DP -an FS -an MQ -an SOR -an MQRankSum -mode SNP --maxGaussians 4 \
--target_titv 2.801 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile primary.recal -tranchesFile primary.recalibrate.tranches -rscriptFile primary.recalibrateSNPs.R

# The most important file is the recalibration report, called recalibrate_SNP.recal, which contains the recalibration data

# now copy all recalibrate* files to your laptop, run the R script, examine the resulting plot and tranches.pdf
# copy files to laptop

# Rscript recalibrateSNPs.R

# explanation for TiTV ratio: https://gatk.broadinstitute.org/hc/en-us/articles/360035531572-Evaluating-the-quality-of-a-germline-short-variant-callset

# Tested the following:
# Lowered the Gausian distributions from 6 to 4 becaue of the low number of SNPs
# Lowered prior from 30 to 10
# Testing all -an QD -an DP -an FS -an MQ -an SOR -an ReadPosRankSum -an MQRankSum and checking distributions
# displaying every 5 tranches from 55 - 100
# Testing with 621 SNPs in training (vqsr2.vcf), 2270 SNPs (vqsr1.vcf) and 2314 (vqsr.vcf)

# The pos SNPs (green) are those which were found in the training sets passed into the VariantRecalibrator step,
# while the neg SNPs (pruple) are those which were found to be furthest away from the learned Gaussians
# and thus given the lowest probability of being true

# Do the annotation dimensions provide a clear separation between the known SNPs (most of which are true)
# and the novel SNPs (most of which are false)?

# ==================== Visualisation of annotations =============================
# This is a guide for visualising annotation in vcf files
# here I have only done it for 6 annotations FS_SOR_MQRS_RPRS_QD_MQ_DP
# for comparision see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
module load bcftools
bcftools query primary_DP3.recode.vcf -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > primary_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt

bgzip -c vqsr1.vcf > vqsr1.vcf.gz
tabix -p vcf vqsr1.vcf.gz
bcftools query vqsr1.vcf.gz -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > vqsr1_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt

scp sries@rackham.uppmax.uu.se:/proj/snic2022-23-656/2bRAD_analysis/fastqs/renamed_bams/primary_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt /Users/xriest/Documents/0_PhD_University_Gothenburg/7_Zostera_marina/Bioinformatics/Seascape_2022/
scp sries@rackham.uppmax.uu.se:/proj/snic2022-23-656/2bRAD_analysis/fastqs/renamed_bams/vqsr1_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt /Users/xriest/Documents/0_PhD_University_Gothenburg/7_Zostera_marina/Bioinformatics/Seascape_2022/


# To plot it run this for all files in R
'#!/usr/bin/env Rscript
library("ggplot2")
library("tidyr")

data <- read.delim("primary_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt", header = F, col.names = c("FS", "SOR", "MQRS", "RPRS", "QD", "MQ", "DP"), na.strings = ".")
data_long <- data %>%                          # Apply pivot_longer function
  pivot_longer(colnames(data)) %>%
  as.data.frame()
ggp3 <- ggplot(data_long, aes(x = value)) +    # Draw histogram & density
  geom_histogram(aes(y = ..density..)) +
  geom_density(col = "#1b98e0", size = 1) +
  facet_wrap(~ name, scales = "free")
ggp3'

'#!/usr/bin/env Rscript
library("ggplot2")
library("tidyr")

data <- read.delim("vqsr1_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt", header = F, col.names = c("FS", "SOR", "MQRS", "RPRS", "QD", "MQ", "DP"), na.strings = ".")
data_long <- data %>%                          # Apply pivot_longer function
  pivot_longer(colnames(data)) %>%
  as.data.frame()
ggp3 <- ggplot(data_long, aes(x = value)) +    # Draw histogram & density
  geom_histogram(aes(y = ..density..)) +
  geom_density(col = "#1b98e0", size = 1) +
  facet_wrap(~ name, scales = "free")
ggp3'

scp sries@rackham.uppmax.uu.se:/proj/snic2022-23-656/2bRAD_analysis/fastqs/renamed_bams/Rplots.pdf /Users/xriest/Documents/0_PhD_University_Gothenburg/7_Zostera_marina/Bioinformatics/Seascape_2022/
# ================================

# applying recalibration:
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
--ts_filter_level 95.0 -mode SNP \
-recalFile primary.recal -tranchesFile primary.recalibrate.tranches -o primary.recal.vcf

# Applying filters
module load bioinfo-tools
module load vcftools
module load picard/2.23.4

vcftools --vcf primary.recal.vcf --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out primary.filt
# After filtering, kept 1053 out of 1053 Individuals
# Outputting VCF file...
# After filtering, kept 7994 out of a possible 63445 Sites

# identifying poorly genotyped individuals
vcftools --vcf primary.filt.recode.vcf --missing-indv --out primary
# look at number of sites genotyped per individual (4th column):
head primary.imiss
# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:

awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' primary.imiss 
# average: 0.096352

cat primary.imiss | awk '$5>0.25' | cut -f 1 > primary.underSequenced
wc -l primary.underSequenced
# 95 under sequenced

grep -F -f primary.underSequenced clonepairs_oldnames.tab | wc -l
# 20

grep -vF -f primary.underSequenced clonepairs_oldnames.tab > primary.clonepairs_good.tab
# applying filter and selecting polymorphic biallelic loci genotyped in 90% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null
# alleles because of mutations in restriction site)


vcftools --vcf primary.filt.recode.vcf --remove primary.underSequenced --max-missing 0.9 --recode-INFO-all --recode --out primary.filt2
# After filtering, kept 959 out of 1053 Individuals
# Outputting VCF file...
# After filtering, kept 6676 out of a possible 7994 Sites


# selecting only polymorphic sites (they all are in denovo pipeline!) and sites with no excess heterozygosity
grep -E "#|0/1|0/0.+1/1|1/1.+0/0" primary.filt2.recode.vcf > primary.polymorphs.vcf
hetfilter.pl vcf=primary.polymorphs.vcf maxhet=0.5 >primary.best.vcf
# 6637 total loci
# 0 dropped because fraction of missing genotypes exceeded 0.5
# 10 dropped because fraction of heterozygotes exceeded 0.5
# 6627 written


# genotypic match between pairs of replicates (the most telling one is the last one, HetsDiscoveryRate - fraction of correctly called heterozygotes; if it is under 90% perhaps use fuzzy genotyping with ANGSD - see above)
repMatchStats.pl vcf=primary.best.vcf replicates=primary.clonepairs_good.tab > primary.hetMatch.txt

#==========================================================
#==========================================================
# SOME CLONEPAIRS WHERE BAD SO HERE I AM REDOING ABOVE VQSR USING ONLY GOOD CLONEPAIRS
#==========================================================
#==========================================================

# many had low HetsDiscoveryRate, so created 3 different filters to test redo replicatesMatch.pl
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.9' | wc -l 
# 73
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.8' | wc -l 
# 93
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.7' | wc -l 
# 106
cat primary.hetMatch.txt | awk -F '\t' '$10>=0.9' | cut -d ':' -f 1  > keep_clonepairs9.tab

grep -F -f keep_clonepairs9.tab clonepairs_oldnames.tab > clonepairs9.tab

# I redid the recalibration with just the good clonepairs (just 73 from 129 were kept)
replicatesMatch.pl vcf=primary_DP3.recode.vcf replicates=clonepairs9.tab hetPairs=1  > vqsr9.vcf
# 63445 total SNPs
# 11551 pass hets and match filters
# 2742 show non-reference alleles
# 2742 have alterantive alleles in at least 1 replicate pair(s)
# 2668 have matching heterozygotes in at least 1 replicate pair(s)
# 2668 polymorphic
# 2668 written

############## I continued with a cutoff of 0.9 => harsh cuttof is recommended

module load bcftools
bgzip -c vqsr9.vcf > vqsr9.vcf.gz
tabix -p vcf vqsr9.vcf.gz
bcftools query vqsr9.vcf.gz -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > vqsr9_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt

vcftools --vcf vqsr9.vcf --TsTv-summary
# Ts/Tv ratio: 2.479

'#!/usr/bin/env Rscript
library("ggplot2")
library("tidyr")

data <- read.delim("vqsr9_FS_SOR_MQRS_RPRS_QD_MQ_DP.txt", header = F, col.names = c("FS", "SOR", "MQRS", "RPRS", "QD", "MQ", "DP"), na.strings = ".")
data_long <- data %>%                          # Apply pivot_longer function
  pivot_longer(colnames(data)) %>%
  as.data.frame()
ggp3 <- ggplot(data_long, aes(x = value)) +    # Draw histogram & density
  geom_histogram(aes(y = ..density..)) +
  geom_density(col = "#1b98e0", size = 1) +
  facet_wrap(~ name, scales = "free")
ggp3'

scp sries@rackham.uppmax.uu.se:/proj/snic2022-23-656/2bRAD_analysis/fastqs/renamed_bams/Rplots.pdf /Users/xriest/Documents/0_PhD_University_Gothenburg/7_Zostera_marina/Bioinformatics/Seascape_2022/

#creating recalibration models
export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"

module load bioinfo-tools
module load picard/2.23.4
module load RStudio
module load R_packages/4.1.1
module load GATK/3.8-0


# OPTION: QD, DP, FS, MQ, SOR, MQRankSum
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar  -T VariantRecalibrator \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
-resource:repmatch,known=true,training=true,truth=true,prior=20 vqsr9.vcf \
-an QD -an DP -an FS -an MQ -an SOR -an MQRankSum -mode SNP --maxGaussians 4 \
--target_titv 2.479 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile second.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrateSNPs.R

# ================================
# copy files to laptop 

# ================================
# applying recalibration:

# I will do this two times
# for the Monitoring file with 99 trench -> continued with this file for the monitoring manuscript

# recalibrate with trench 99 - MONITORING
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
--ts_filter_level 99.0 -mode SNP \
-recalFile second.recal -tranchesFile recalibrate.tranches -o recal.vcf

# recalibrate with trench 95 -  SEASCAPE
java -Djava.io.tmpdir=$HOME -jar $GATK_HOME/GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $REFERENCE_GENOME -input primary_DP3.recode.vcf -nt 1 \
--ts_filter_level 95.0 -mode SNP \
-recalFile second.recal -tranchesFile recalibrate.tranches -o recal.vcf

# Applying filters
module load bioinfo-tools
module load vcftools

# applying filter and selecting polymorphic biallelic loci genotyped in 90% or more individuals
# (harsh genotyping rate cutoff is strongly recommended for best quality and to avoid RAD loci affected by null
# alleles because of mutations in restriction site)

vcftools --vcf recal.vcf --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --out filt
# After filtering, kept 1077 out of 1077 Individuals
# Outputting VCF file...
# After filtering, kept 13243 out of a possible 64345 Sites


# identifying poorly genotyped individuals
vcftools --vcf filt.recode.vcf --missing-indv --out out2
# look at number of sites genotyped per individual (4th column):
head out2.imiss

# see if some samples are much lower in the number of sites than others
# for example, if you want to remove samples showing less than 40000 sites:
awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' out2.imiss 
## 99 - MONITORING 
# 0.147764


cat out2.imiss | awk '$5>0.42' | cut -f 1 > underSequenced2
wc -l underSequenced2
# 69 underSequenced2

grep -F -f underSequenced2 clonepairs_renamed.tab | wc -l
# 13 underSequenced2

grep -vF -f underSequenced2 clonepairs_renamed.tab > clonepairs_good2.tab

vcftools --vcf filt.recode.vcf --remove underSequenced2 --max-missing 0.9 --recode-INFO-all --recode --out filt2
# After filtering, kept 1009 out of 1077 Individuals
# Outputting VCF file...
# After filtering, kept 9069 out of a possible 13243 Site


grep -n '#C' filt2.recode.vcf | cut -f10-1000 | sed 's/\t/\n/g' > indv_filt2_vcf


# selecting only polymorphic sites (they all are in denovo pipeline!) and sites with no excess heterozygosity
grep -E "#|0/1|0/0.+1/1|1/1.+0/0" filt2.recode.vcf >polymorphs.vcf
hetfilter.pl vcf=polymorphs.vcf maxhet=0.5 >best.vcf

# 9041 total loci
# 0 dropped because fraction of missing genotypes exceeded 0.5
# 35 dropped because fraction of heterozygotes exceeded 0.5
# 9006 written

grep -n '#C' best.vcf | cut -f10-1000 | sed 's/\t/\n/g' > indv_best_vcf

#---------------
# Final touches

#-------------------
# KEEPING replicates
#-----------------------
thinner.pl infile=best.vcf criterion=maxAF >thinMaxaf.vcf
# 9006 total loci
# 5771 loci selected

# genotypic match between pairs of replicates (the most telling one is the last one, HetsDiscoveryRate - fraction of correctly called heterozygotes; if it is under 90% perhaps use fuzzy genotyping with ANGSD - see above)
repMatchStats.pl vcf=thinMaxaf.vcf replicates=clonepairs_good2.tab > hetMatch.txt
awk -F '\t' '$10<0.9' hetMatch.txt | wc -l
# 70

#-------------------------------------------------
### little bit further I had problems with the renamed files -> therfore I kept the old labelling for now with "oldnames"
vcftools --vcf thinMaxaf.vcf --missing-indv --out minDP
# After filtering, kept 1009 out of 1009 Individuals
# Outputting Individual Missingness
# After filtering, kept 5771 out of a possible 5771 Sites

awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' minDP.imiss
# 0.0243533

#----------
# 25% MISSINGNES
# -----------
cat minDP.imiss | awk '$5>0.25' | cut -f 1  > ind_miss_25
wc -l ind_miss_25
# 8 ind_miss_25


vcftools --vcf thinMaxaf.vcf --remove ind_miss_25 --max-missing 0.9 --recode-INFO-all --recode --out all_miss25
# After filtering, kept 1002 out of 1009 Individuals
# Outputting VCF file...
# After filtering, kept 5771 out of a possible 5771 Sites

vcftools --vcf thinMaxaf.vcf --remove ind_miss_25 --max-missing 0.9 --maf 0.01 --recode-INFO-all --recode --out all_miss25_maf
# After filtering, kept 1002 out of 1009 Individuals
# Outputting VCF file...
# After filtering, kept 581 out of a possible 5771 Sites

#----------
# 5% MISSINGNES
# -----------
cat minDP.imiss | awk '$5>0.05' | cut -f 1  > ind_miss_05
wc -l ind_miss_05
# 153 ind_miss_05

vcftools --vcf thinMaxaf.vcf --remove ind_miss_05 --max-missing 0.9 --recode-INFO-all --recode --out all_miss05
# After filtering, kept 857 out of 1009 Individuals
# Outputting VCF file...
# After filtering, kept 5771 out of a possible 5771 Sites


vcftools --vcf thinMaxaf.vcf --remove ind_miss_05 --max-missing 0.9 --maf 0.01 --recode-INFO-all --recode --out all_miss05_maf
# After filtering, kept 857 out of 1009 Individuals
# Outputting VCF file...
# After filtering, kept 575 out of a possible 5771 Sites


# het match between replicates
grep -F -f indv clonepairs_good2.tab | grep -vF -f ind_miss_25 > clonepairs_renamed_25.tab
grep -F -f indv clonepairs_good2.tab | grep -vF -f ind_miss_05 > clonepairs_renamed_05.tab

repMatchStats.pl vcf=all_miss25.recode.vcf replicates=clonepairs_renamed_25.tab > hetMatch_all_miss25.txt
repMatchStats.pl vcf=all_miss05.recode.vcf replicates=clonepairs_renamed_05.tab > hetMatch_all_miss05.txt
awk -F '\t' '$10<0.9' hetMatch.txt | wc -l
# 70 -7 = 63


awk -F '\t' '$10<0.9' hetMatch_all_miss25.txt | wc -l
# 69 -7 = 62

awk -F '\t' '$10<0.9' hetMatch_all_miss05.txt | wc -l
# 40 -7 = 33

##########################################

# SPLIT AND SORT
#-------------------------------------------------
# renaming individuals
#-------------------------------------------------
# rename now the vcf files!
# RENAMING all_miss vcf file

#-------------------------------------------------
# all_miss25
#-------------------------------------------------
grep -n '#C' all_miss25.recode.vcf | cut -f10-1500 | sed 's/\t/\n/g' > indv_all_miss25
sed -i.backup 's/YST-19/STE-18-rep/g; s/BJO-4-Rep/ALA-05-rep/g; s/FUR-20-Rep/HOR-11-rep/g; s/GRO-8-Rep/HOG-12-rep/g; s/HOG-12-Rep/GRO-08-rep/g; s/HOR-11-Rep/FUR-20-rep/g;
s/ALA-5-Rep/BJO-04-rep/g; s/STE-18-Rep/YST-19/g; s/VF-3360-BUR-11/VF-3360-BUR-11-rep/g; s/VJ-3470-LOM-20-rep/VJ-3470-LOM-20/g; sx-\([1-9]_\)x-0\1xg; sx-\([1-9]-\)x-0\1xg; s/_S.*$//g; s/a$//g; s/_R1//g; s/KO\-DON/KOD/g; s/repl/rep/g; s/Rep/rep/g; s/rep1/rep/g; s/rep2/rep/g; s/NYN-12-M/NYN-12/g; s/NYN-13-M/NYN-13/g; 
s/NYN-14-M/NYN-14/g; s/NYN-14-R/NYN-14-rep/g; s/FLA03/FLA-03/g; s/FLA10/FLA-10/g; s/KOC15/KOC-15/g; s/KOC8/KOC-08/g; s/NYC07/NYC-07/g; s/NYC10/NYC-10/g; 
s/RAM04/RAM-04/g; s/RAM06/RAM-06/g; s/STY07/STY-07/g; s/STY12/STY-12/g; s/TAN11/TAN-11/g; s/TAN15/TAN-15/g; s/TJB08/TJB-08/g; s/TJB19/TJB-19/g; s/VAT13/VAT-13/g; s/VAT15/VAT-15/g; s/MS22/M22/g; s/ANG/YU-0000-ANG/g' indv_all_miss25

paste indv_all_miss25.backup indv_all_miss25 > rename_indv_all_miss25

module load bcftools
bcftools reheader -s rename_indv_all_miss25 all_miss25.recode.vcf > all_miss25_renamed.vcf

# split datasets
grep -E "FLA|GAS|KOC|KOD|M22|MS22|NYC|RAM|STY|TAN|TJB|VAT" indv_all_miss25 | sort > maru_indv_all_25
grep -v "FLA\|KOC\|M22\|MS22\|NYC\|RAM\|STY\|TAN\|TJB\|VAT" indv_all_miss25 | sort > steffi_indv_all_25
grep -E "KOD|STE|GOT|GRO|HOG|ALA|YST|FUR|KAR|HOR|KRA|KLI|SLI|NYN|BJO" indv_all_miss25 | sort > steffi_monitoring_indv_all_25
# ind list without replicates
grep 'rep\|repl\|Rep' steffi_indv_all_25 > steffi_indv_uniq_25

#-------------------------------------------------
# all_miss25_maf
#-------------------------------------------------
grep -n '#C' all_miss25_maf.recode.vcf | cut -f10-1500 | sed 's/\t/\n/g' > indv_all_miss25_maf
sed -i.backup 's/YST-19/STE-18-rep/g; s/BJO-4-Rep/ALA-05-rep/g; s/FUR-20-Rep/HOR-11-rep/g; s/GRO-8-Rep/HOG-12-rep/g; s/HOG-12-Rep/GRO-08-rep/g; s/HOR-11-Rep/FUR-20-rep/g;
s/ALA-5-Rep/BJO-04-rep/g; s/STE-18-Rep/YST-19/g; s/VF-3360-BUR-11/VF-3360-BUR-11-rep/g; s/VJ-3470-LOM-20-rep/VJ-3470-LOM-20/g; sx-\([1-9]_\)x-0\1xg; sx-\([1-9]-\)x-0\1xg; s/_S.*$//g; s/a$//g; s/_R1//g; s/KO\-DON/KOD/g; s/repl/rep/g; s/Rep/rep/g; s/rep1/rep/g; s/rep2/rep/g; s/NYN-12-M/NYN-12/g; s/NYN-13-M/NYN-13/g; 
s/NYN-14-M/NYN-14/g; s/NYN-14-R/NYN-14-rep/g; s/FLA03/FLA-03/g; s/FLA10/FLA-10/g; s/KOC15/KOC-15/g; s/KOC8/KOC-08/g; s/NYC07/NYC-07/g; s/NYC10/NYC-10/g; 
s/RAM04/RAM-04/g; s/RAM06/RAM-06/g; s/STY07/STY-07/g; s/STY12/STY-12/g; s/TAN11/TAN-11/g; s/TAN15/TAN-15/g; s/TJB08/TJB-08/g; s/TJB19/TJB-19/g; s/VAT13/VAT-13/g; s/VAT15/VAT-15/g; s/MS22/M22/g; s/ANG/YU-0000-ANG/g' indv_all_miss25_maf

paste indv_all_miss25_maf.backup indv_all_miss25_maf > rename_indv_all_miss25_maf

bcftools reheader -s rename_indv_all_miss25_maf all_miss25_maf.recode.vcf > all_miss25_maf_renamed.vcf

# split datsets
grep -E "FLA|GAS|KOC|KOD|M22|MS22|NYC|RAM|STY|TAN|TJB|VAT" indv_all_miss25_maf | sort > maru_indv_all_25_maf
grep -v "FLA\|KOC\|M22\|MS22\|NYC\|RAM\|STY\|TAN\|TJB\|VAT" indv_all_miss25_maf | sort > steffi_indv_all_25_maf
grep -E "KOD|STE|GOT|GRO|HOG|ALA|YST|FUR|KAR|HOR|KRA|KLI|SLI|NYN|BJO" indv_all_miss25_maf | sort > steffi_monitoring_indv_all_25_maf
# ind list without replicates
grep 'rep\|repl\|Rep' steffi_indv_all_25_maf > steffi_indv_uniq_25_maf

#-------------------------------------------------
# all_miss05
#-------------------------------------------------
grep -n '#C' all_miss05.recode.vcf | cut -f10-1500 | sed 's/\t/\n/g' > indv_all_miss05
sed -i.backup 's/YST-19/STE-18-rep/g; s/BJO-4-Rep/ALA-05-rep/g; s/FUR-20-Rep/HOR-11-rep/g; s/GRO-8-Rep/HOG-12-rep/g; s/HOG-12-Rep/GRO-08-rep/g; s/HOR-11-Rep/FUR-20-rep/g;
s/ALA-5-Rep/BJO-04-rep/g; s/STE-18-Rep/YST-19/g; s/VF-3360-BUR-11/VF-3360-BUR-11-rep/g; s/VJ-3470-LOM-20-rep/VJ-3470-LOM-20/g; sx-\([1-9]_\)x-0\1xg; sx-\([1-9]-\)x-0\1xg; s/_S.*$//g; s/a$//g; s/_R1//g; s/KO\-DON/KOD/g; s/repl/rep/g; s/Rep/rep/g; s/rep1/rep/g; s/rep2/rep/g; s/NYN-12-M/NYN-12/g; s/NYN-13-M/NYN-13/g; 
s/NYN-14-M/NYN-14/g; s/NYN-14-R/NYN-14-rep/g; s/FLA03/FLA-03/g; s/FLA10/FLA-10/g; s/KOC15/KOC-15/g; s/KOC8/KOC-08/g; s/NYC07/NYC-07/g; s/NYC10/NYC-10/g; 
s/RAM04/RAM-04/g; s/RAM06/RAM-06/g; s/STY07/STY-07/g; s/STY12/STY-12/g; s/TAN11/TAN-11/g; s/TAN15/TAN-15/g; s/TJB08/TJB-08/g; s/TJB19/TJB-19/g; s/VAT13/VAT-13/g; s/VAT15/VAT-15/g; s/MS22/M22/g; s/ANG/YU-0000-ANG/g' indv_all_miss05

paste indv_all_miss05.backup indv_all_miss05 > rename_indv_all_miss05

module load bcftools
bcftools reheader -s rename_indv_all_miss05 all_miss05.recode.vcf > all_miss05_renamed.vcf

grep -E "FLA|GAS|KOC|KOD|M22|MS22|NYC|RAM|STY|TAN|TJB|VAT" indv_all_miss05 | sort > maru_indv_all_05
grep -v "FLA\|KOC\|M22\|MS22\|NYC\|RAM\|STY\|TAN\|TJB\|VAT" indv_all_miss25 | sort > steffi_indv_all_05
grep -E "KOD|STE|GOT|GRO|HOG|ALA|YST|FUR|KAR|HOR|KRA|KLI|SLI|NYN|BJO" indv_all_miss05 | sort > steffi_monitoring_indv_all_05

# ind list without replicates
grep 'rep\|repl\|Rep' steffi_indv_all_05 > steffi_indv_uniq_05

#-------------------------------------------------
# all_miss05_maf
#-------------------------------------------------
grep -n '#C' all_miss05_maf.recode.vcf | cut -f10-1500 | sed 's/\t/\n/g' > indv_all_miss05_maf
sed -i.backup 's/YST-19/STE-18-rep/g; s/BJO-4-Rep/ALA-05-rep/g; s/FUR-20-Rep/HOR-11-rep/g; s/GRO-8-Rep/HOG-12-rep/g; s/HOG-12-Rep/GRO-08-rep/g; s/HOR-11-Rep/FUR-20-rep/g;
s/ALA-5-Rep/BJO-04-rep/g; s/STE-18-Rep/YST-19/g; s/VF-3360-BUR-11/VF-3360-BUR-11-rep/g; s/VJ-3470-LOM-20-rep/VJ-3470-LOM-20/g; sx-\([1-9]_\)x-0\1xg; sx-\([1-9]-\)x-0\1xg; s/_S.*$//g; s/a$//g; s/_R1//g; s/KO\-DON/KOD/g; s/repl/rep/g; s/Rep/rep/g; s/rep1/rep/g; s/rep2/rep/g; s/NYN-12-M/NYN-12/g; s/NYN-13-M/NYN-13/g; 
s/NYN-14-M/NYN-14/g; s/NYN-14-R/NYN-14-rep/g; s/FLA03/FLA-03/g; s/FLA10/FLA-10/g; s/KOC15/KOC-15/g; s/KOC8/KOC-08/g; s/NYC07/NYC-07/g; s/NYC10/NYC-10/g; 
s/RAM04/RAM-04/g; s/RAM06/RAM-06/g; s/STY07/STY-07/g; s/STY12/STY-12/g; s/TAN11/TAN-11/g; s/TAN15/TAN-15/g; s/TJB08/TJB-08/g; s/TJB19/TJB-19/g; s/VAT13/VAT-13/g; 
s/VAT15/VAT-15/g; s/MS22/M22/g; s/ANG/YU-0000-ANG/g' indv_all_miss05_maf

paste indv_all_miss05_maf.backup indv_all_miss05_maf > rename_indv_all_miss05_maf

module load bcftools
bcftools reheader -s rename_indv_all_miss05_maf all_miss05_maf.recode.vcf > all_miss05_maf_renamed.vcf

grep -E "FLA|GAS|KOC|KOD|M22|MS22|NYC|RAM|STY|TAN|TJB|VAT" indv_all_miss05_maf | sort > maru_indv_all_05_maf
grep -v "FLA\|KOC\|M22\|MS22\|NYC\|RAM\|STY\|TAN\|TJB\|VAT" indv_all_miss25_maf | sort > steffi_indv_all_05_maf
grep -E "KOD|STE|GOT|GRO|HOG|ALA|YST|FUR|KAR|HOR|KRA|KLI|SLI|NYN|BJO" indv_all_miss05_maf | sort > steffi_monitoring_indv_all_05_maf

# ind list without replicates
grep 'rep\|repl\|Rep' steffi_indv_all_05_maf > steffi_indv_uniq_05_maf

###### VCF FILES ######

######################
###### MISS 025 ######
#----------------------------
vcftools --vcf all_miss25_renamed.vcf --keep steffi_monitoring_indv_all_25 --recode --recode-INFO-all --out zostera_monitoring_230504_steffi_miss25_99 # Steffis
# After filtering, kept 260 out of 1002 Individuals
# Outputting VCF file...
# After filtering, kept 5771 out of a possible 5771 Sites

######################
###### MISS 05 ######
#----------------------------
vcftools --vcf all_miss05_renamed.vcf --keep steffi_monitoring_indv_all_05 --recode --recode-INFO-all --out zostera_monitoring_230504_steffi_miss05_99 # Steffis
# After filtering, kept 203 out of 857 Individuals
# Outputting VCF file...
# After filtering, kept 5771 out of a possible 5771 Sites

################################################
###### MISS MAF 25 and 05 for all my indv ######
#----------------------------
vcftools --vcf all_miss25_maf_renamed.vcf --keep steffi_indv_all_25_maf --recode --recode-INFO-all --out zostera_230504_steffi_miss25_maf_99 # Steffis
# After filtering, kept 774 out of 1002 Individuals
# Outputting VCF file...
# After filtering, kept 581 out of a possible 581 Sites


vcftools --vcf all_miss05_maf_renamed.vcf --keep steffi_indv_all_05_maf --recode --recode-INFO-all --out zostera_230504_steffi_miss05_maf_99 # Steffis
# After filtering, kept 653 out of 857 Individuals
# Outputting VCF file...
# After filtering, kept 575 out of a possible 575 Sites


####################################
###### Unique data sets #####
vcftools --vcf all_miss25_renamed.vcf --remove steffi_indv_uniq_25 --recode --recode-INFO-all --out zostera_230414_230504_miss25_uniq_99 # uniq
# After filtering, kept 900 out of 1002 Individuals
# Outputting VCF file...
# After filtering, kept 5771 out of a possible 5771 Sites

vcftools --vcf all_miss25_maf_renamed.vcf --remove steffi_indv_uniq_25_maf --recode --recode-INFO-all --out zostera_230504_steffi_miss25_maf_uniq_99
# After filtering, kept 900 out of 1002 Individuals
# Outputting VCF file...
# After filtering, kept 581 out of a possible 581 Sites

vcftools --vcf all_miss05_renamed.vcf --remove steffi_indv_uniq_05 --recode --recode-INFO-all --out zostera_230504_steffi_miss05_uniq_99 # uniq
# After filtering, kept 769 out of 857 Individuals
# Outputting VCF file...
# After filtering, kept 5771 out of a possible 5771 Sites

vcftools --vcf all_miss05_maf_renamed.vcf --remove steffi_indv_uniq_05_maf --recode --recode-INFO-all --out zostera_230504_steffi_miss05_maf_uniq_99 # uniq
# After filtering, kept 769 out of 857 Individuals
# Outputting VCF file...
# After filtering, kept 575 out of a possible 575 Sites

#----------------------------------------------
#SORT
# sort individuals according to a list
#----------------------------------------------

# pop order
echo "KIE
KOD
GAS
KAV
VHA
STE
NOV
SBR
GOT
GRO
HOG
ALA
BAR
LOM
KUR
YST
SKI
AHU
FUR
BJN
BLA
KAL
KAR
OSK
BUR
BAD
LJU
KLI
HOR
SLI
EKO
KRA
NYN
GAL
BJO
RAU
SAR
HAN
ING
TAL" > pop_order

# reorders the file

#-------- for 25 monitoring
:> steffi_monitoring_indv_all_sorted_25
while read TESTLINE; do
  grep "${TESTLINE}" steffi_monitoring_indv_all_25 >> steffi_monitoring_indv_all_sorted_25
done < pop_order

#-------- for 05 monitoring
:> steffi_monitoring_indv_all_sorted_05
while read TESTLINE; do
  grep "${TESTLINE}" steffi_monitoring_indv_all_05 >> steffi_monitoring_indv_all_sorted_05
done < pop_order

module load bioinfo-tools
module load bcftools
module load vcftools

#### ALL INDIVIDUALS VCF FILES
bcftools view -S steffi_indv_all_sorted_25 --force-samples zostera_230504_steffi_miss25_99.recode.vcf -O v -o zostera_230504_steffi_miss25_sorted_99.vcf
bcftools view -S steffi_indv_all_sorted_05 --force-samples zostera_230504_steffi_miss05_99.recode.vcf -O v -o zostera_230504_steffi_miss05_sorted_99.vcf
bcftools view -S steffi_indv_all_sorted_25_maf --force-samples zostera_230504_steffi_miss25_maf_99.recode.vcf -O v -o zostera_230504_steffi_miss25_maf_sorted_99.vcf
bcftools view -S steffi_indv_all_sorted_05_maf --force-samples zostera_230504_steffi_miss05_maf_99.recode.vcf -O v -o zostera_230504_steffi_miss05_maf_sorted_99.vcf

#### MONITORING VCF FILES
bcftools view -S steffi_monitoring_indv_all_sorted_25 --force-samples zostera_monitoring_230504_steffi_miss25_99.recode.vcf -O v -o zostera_monitoring_230504_steffi_miss25_sorted_99.vcf
bcftools view -S steffi_monitoring_indv_all_sorted_05 --force-samples zostera_monitoring_230504_steffi_miss05_99.recode.vcf -O v -o zostera_monitoring_230504_steffi_miss05_sorted_99.vcf

###############
###   END   ###
###############

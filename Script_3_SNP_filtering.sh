#####################################
#
# Part 3: SNP filtering (VCFtools)
# January 2025
#
# Code Contributor: Cassandra E. & Yue Y
# Biodiversity Research Center, UBC
#
#####################################



# This code includes the following:
# 3.1: SNP filtering (VCFtools)
#      quality, minor allele counts, % found accross individuals
#      filter out individuals with lots of missing data
#      A script to count the number of potential genotyping errors due to low read depth




#####################################
#  3.1 SNP filtering (VCFtools)
#####################################

# -- raw VCF file with 40 samples (10 individuals/populations)
# Note: Four populations in this VCF are: GEN, BARD, KL, PC

cd /scratch/celphin/GBS_workshop/0_Raw_data


# -- load vcftools/bcftools
module load StdEnv/2020 
module load vcftools/0.1.16
module load gnuplot/5.4.2
module load gcc/9.3.0
module load bcftools/1.16
module load plink/1.9b_6.21-x86_64


# -- Count number of sites (SNPs + others)
grep -v '^#' TotalRawSNPs_subset_40samples.vcf | wc -l
# 559,815 Variants



# ---------------------------------------------------
#  filtering SNP step 1 : start new computing session 
# ---------------------------------------------------
tmux new-session -s test
tmux attach-session -t test

salloc -c1 --time 03:00:00 --mem 30G --account def-rieseber #(change this to your account)



cd /scratch/celphin/GBS_workshop/3_SNP_filtering

# Copy the raw VCF file to the current directory
cp /scratch/celphin/GBS_workshop/0_Raw_data/TotalRawSNPs_subset_40samples.vcf .

ls



# ---------------------------------------------------
#  filtering SNP step 2 : HWE filter (high heterzygosity)
# ---------------------------------------------------
# convert the vcf file to Plink
vcftools --vcf TotalRawSNPs_subset_40samples.vcf --plink --out TotalRawSNPs_plink

# find SNPs heterozygous in all or 90% of individuals 
plink --file TotalRawSNPs_plink --hardy --out TotalRawSNPs_hardy

more TotalRawSNPs_hardy.hwe

# find SNPs with high heterozygosity
mawk '$7 > 0.4' TotalRawSNPs_hardy.hwe | cut -f1 > HighHetSNPs_40.txt 
mawk '$7 > 0.5' TotalRawSNPs_hardy.hwe | cut -f1 > HighHetSNPs_50.txt 
mawk '$7 > 0.6' TotalRawSNPs_hardy.hwe | cut -f1 > HighHetSNPs_60.txt 
mawk '$7 > 0.7' TotalRawSNPs_hardy.hwe | cut -f1 > HighHetSNPs_70.txt 
mawk '$7 > 0.8' TotalRawSNPs_hardy.hwe | cut -f1 > HighHetSNPs_80.txt
mawk '$7 > 0.9' TotalRawSNPs_hardy.hwe | cut -f1 > HighHetSNPs_90.txt

wc -l HighHetSNPs_40.txt
#8017 HighHetSNPs_40.txt
wc -l HighHetSNPs_50.txt
#5074 HighHetSNPs_50.txt
wc -l HighHetSNPs_60.txt
#3813 HighHetSNPs_60.txt
wc -l HighHetSNPs_70.txt
#2969 HighHetSNPs_70.txt
wc -l HighHetSNPs_80.txt
#2342 HighHetSNPs_80.txt
wc -l HighHetSNPs_90.txt
#2000 HighHetSNPs_90.txt


# get list of SNP IDs
awk '{print $2}' HighHetSNPs_60.txt > HighHetSNPs_60_list
head HighHetSNPs_60_list

# Change names into 2 columns
sed 's/:/\t/g' HighHetSNPs_60_list > HighHetSNPs_60_list.txt
head HighHetSNPs_60_list.txt


# filter out this list of SNPs that are 60% het across all populations
VCF="TotalRawSNPs_subset_40samples.vcf"
LIST="HighHetSNPs_60_list.txt"

vcftools --vcf $VCF --exclude-positions $LIST --recode --recode-INFO-all --out TotalRawSNPs_rmhet
# After filtering, kept 556,003 out of a possible 559,815 Sites




# ---------------------------------------------------
#  filtering SNP step 3 : Rename samples in vcf file
# ---------------------------------------------------
# -- https://www.biostars.org/p/279195/ 

bcftools query -l TotalRawSNPs_rmhet.recode.vcf > sample_names.txt

#sed 's/,/ /g' M_caridnalis_samples_renamed.csv > M_caridnalis_samples_renamed.txt
#bcftools reheader -s M_caridnalis_samples_renamed.txt -o Mimulus_timeseries_filtered_variants_rename.vcf Mimulus_timeseries_filtered_variants.vcf
#bcftools query -l Mimulus_timeseries_filtered_variants_rename.vcf



# ---------------------------------------------------
#  filtering SNP step 4 : Site-level Filtering
# ---------------------------------------------------
# -- ON the following:

# -> high quality
# -> remove indels
# -> biallelic
# -> missing data < 90%
# -> monomorphic
# -> LD
# -> minor allele freq of 0.01 (0.01 * 40 * 2 = 8 alleles)

# -- Code
vcftools --vcf TotalRawSNPs_rmhet.recode.vcf \
--minGQ 30 \
--remove-indels \
--max-alleles 2 \
--max-missing 0.9 \
--min-alleles 2 \
--thin 300 \ 
# filteirng for LD (1 SNP to every 300bp)
--maf 0.01 \
--recode \
--recode-INFO-all \
--out Cassiope_noMERhb

# After filtering, kept 10,857 out of a possible 556,003 Sites



# ---------------------------------------------------
#  filtering SNP step 5 : Sample-level Filtering
# ---------------------------------------------------

# -- list the amount of missing data per individual - find indivdiuals with no reads mapped
vcftools --vcf Cassiope_noMERhb.recode.vcf --missing-indv --out Cassiope_noMER


# -- filter out the individuals with greater than 10% missing SNPs 
mawk '$5 > 0.1' Cassiope_noMER.imiss | cut -f1 > lowDP.indv
vcftools --vcf Cassiope_noMERhb.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out Cassiope_noMER_r10i
#After filtering, kept 39 out of 40 Individuals
#After filtering, kept 10857 out of a possible 10857 Sites



# END



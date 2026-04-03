#####################################
#
# Part 4: SNP filtering (VCFtools)
# January 2025
#
# Code Contributor: Cassandra E. & Yue Y
# Biodiversity Research Center, UBC
#
#####################################

# This code includes the following:

# -- OVERVIEW GBS SNP FILTERING STEPS

# -- Step 1: Check your data ( no code )
# -- Step 2: Extract SNPs/INDELs
# -- Step 3: Plot out parameters for INFO field 
# -- Step 4: GATK hard filter on INFO field for reliable SNPs (NOT related to N sample size) AND general loose ExcessHet filter
# -- Step 5: Check per sample DP and GQ before sample-level filters
# -- Step 6: GATK filter on SAMPLE field for SAMPLE-DP and SAMPLE-GQ (NOT related to N sample size)
# -- Step 7: Remove sample with low coverage
# -- Step 8: Recode INFO field
# -- Step 9: GATK filter on SAMPLE field pop specific parameters (YES related to N sample size), eg. MAF/MAC/missing rate etc...





#--------------------------
#  Step 2: Extract SNPs
#--------------------------

nano Filter_Step1.sh
# ----------- Filter_Step1.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac 
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=2:00:00

cd /home/yueyu/scratch/GBS/VCF_DEB_FILTERED

module load StdEnv/2023
module load gatk/4.6.1.0

unset JAVA_TOOL_OPTIONS

#SNP
gatk --java-options "-Xmx8g" SelectVariants \
  -V /home/yueyu/scratch/GBS/VCF_DEB/Raw_VCF.vcf.gz \
  --select-type-to-include SNP \
  --exclude-filtered \
  -O SNP_FILTERED.vcf.gz
# ----------- Filter_Step1.sh (END)---------




#--------------------------
#  Step 3: Plot SNP quality (DO NOT RUN)
#--------------------------
module load StdEnv/2023  
module load gcc/12.3
module load bcftools/1.19


#Extract info and plot in R 
bcftools query SNP_v20250415.vcf.gz -f '%ExcessHet\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > SNP_v20250415_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt
wc -l SNP_v20250415_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt
# 2,342,801 (2 million SNPs - DEB)


module load StdEnv/2023
module load r/4.4.0

R
# ------------- in R -------------
library(data.table)

# Load DEB SNP table
snps <- fread("SNP_v20250415_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt", na.strings = c(".", "NA"), colClasses = "numeric", header = F)

# Load PRA SNP table
# snps <- fread("PRA_SNP_v20250415_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt", na.strings = c(".", "NA"), colClasses = "numeric", header = F)

str(snps)
colnames(snps) <- c("ExcessHet","FS","SOR","MQRankSum","ReadPosRankSum","QD","MQ","DP")

# !!! IMPORTANT!! change saved PNG names

# -- ExcessHet
png(filename = "ExcessHet_v20150415.png")
dExcessHet <- density(snps$ExcessHet,na.rm=T)
plot(dExcessHet,main="ExcessHet distribution", xlab="ExcessHet")
dev.off()

# -- FS
# or can plot FS with a log base 10 scale
png(filename = "FS_v20150415.png")
dFS <- density(snps$FS,na.rm=T)
plot(dFS,main="FS distribution", xlab="FS")
dev.off()

# -- SOR
png(filename = "SOR_v20150415.png")
dSOR <- density(snps$SOR,na.rm=T)
plot(dSOR,main="SOR distribution", xlab="SOR")
dev.off()

# -- MQRankSum
png(filename = "MQRankSum_v20150415.png")
dMQRankSum <- density(snps$MQRankSum,na.rm=T)
plot(dMQRankSum,main="MQRankSum distribution", xlab="MQRankSum")
dev.off()

# -- ReadPosRankSum
png(filename = "ReadPosRankSum_v20150415.png")
dReadPosRankSum <- density(snps$ReadPosRankSum,na.rm=T)
plot(dReadPosRankSum,main="ReadPosRankSum distribution", xlab="ReadPosRankSum")
dev.off()

# -- QD
png(filename = "QD_v20150415.png")
dQD <- density(snps$QD,na.rm=T)
plot(dQD,main="QD distribution", xlab="QD")
dev.off()

# -- MQ
png(filename = "MQ_v20150415.png")
dMQ <- density(snps$MQ,na.rm=T)
plot(dMQ,main="MQ distribution", xlab="MQ")
dev.off()

# -- DP
png(filename = "DP_v20150415.png")
dDP <- density(snps$DP,na.rm=T)
plot(dDP,main="DP distribution", xlab="DP")
dev.off()



#--------------------------
#  Step 4: Hard filter
#--------------------------

nano Filter_Step2.sh
# ----------- Filter_Step2.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1:00:00

cd /home/yueyu/scratch/GBS/VCF_DEB_FILTERED

module load StdEnv/2023
module load gatk/4.6.1.0

unset JAVA_TOOL_OPTIONS

gatk --java-options "-Xmx8g" VariantFiltration \
    -V SNP_FILTERED.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSumNeg12.5" \
    -filter "ExcessHet > 54.69" --filter-name "ExcessHet5469" \
    -O SNP_INFO_MARKED.vcf.gz


# ADD CODE TO EXCLUDE THE FILTERED
gatk --java-options "-Xmx8g" SelectVariants \
  -V SNP_INFO_MARKED.vcf.gz \
  --exclude-filtered \
  -O SNP_INFO_FILTERED.vcf.gz

# ----------- Filter_Step2.sh (END)---------


# --------------------
# Step 5: Check per sample DP and GQ before sample-level filters
# --------------------

# Example
# ===========
#  (Genotype) GQ, mean ~ 58 
# ===========
# ===========
#  (Genotype) DP, mean ~ 22
# ===========


# --------------------
# Step 6: BCFTOOLS on SAMPLE field for SAMPLE-DP and SAMPLE-GQ (Zhe code, very good!)
# --------------------



nano Filter_Step3.sh
# ----------- Filter_Step3.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00

cd /home/yueyu/scratch/GBS/VCF_DEB_FILTERED

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6

bcftools filter \
    -e '(GQ<30 | FORMAT/DP<10)' \
    SNP_INFO_FILTERED.vcf.gz \
    -S . \
    -s "PASS" \
    -Oz -o INFO_GENO_FILTERED_single_sep.vcf.gz

tabix -p vcf INFO_GENO_FILTERED_single_sep.vcf.gz

# ----------- Filter_Step3.sh (END)---------


# --------------------
# Step 7: FILTER MAX MISSINGNESS (site level)
# --------------------
nano Filter_Step4.sh
# ----------- Filter_Step4.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00

cd /home/yueyu/scratch/GBS/VCF_DEB_FILTERED

module load StdEnv/2023
module load gatk/4.6.1.0

unset JAVA_TOOL_OPTIONS

gatk --java-options "-Xmx18g" SelectVariants \
            -V SNP_INFO_GENO_FILTERED_single_sep.vcf.gz \
            --remove-unused-alternates \
            --restrict-alleles-to BIALLELIC \
            --max-nocall-fraction 0.3 \
            -O SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz

# ----------- Filter_Step4.sh (END)---------


# --------------------
# Step 8: REMOVE SAMPLES THAT ARE OF LOW COVERAGE -- YES, do it after filter for missing rate (of the sites)
# -------------------- 

module load StdEnv/2020 
module load vcftools/0.1.16

cd /home/yueyu/scratch/GBS/VCF_DEB_FILTERED
vcftools --gzvcf SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz --missing-indv --out sample_missing_ratio
awk '$5 > 0.3' sample_missing_ratio.imiss



# ------------------------------------
# Step 9: RECODE INFO + AF filter
# -----------------------------------


nano Filter_Step5.sh
# ----------- Filter_Step5.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=2:00:00

cd /home/yueyu/scratch/GBS/VCF_DEB_FILTERED

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6

# not needed and recode information (INFO - AN, AC, AF, DP need to recode for downstream)
# Remove other INFO field that are not accurate anymore and can not be re-calculated after filtering (eg. require initial seq info)
# ^ means to keep these INFO fields but drop the rest

bcftools annotate -x "^INFO/AN,INFO/AC,INFO/AF,INFO/DP" SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz -Ou | bcftools +fill-tags  -O z -o SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED.vcf.gz -- -t AN,AC,AF
tabix -p vcf SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED.vcf.gz

# -- AF filter
bcftools view \
  -i 'AF>=0.03 && AF<=0.97' \
  -O z \
  -o SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED_AF003.vcf.gz \
  SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED.vcf.gz

tabix -p vcf SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED_AF003.vcf.gz

# ----------- Filter_Step5.sh (END) ---------

# Check number filtered SNPs left
bcftools view -H SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED_AF003.vcf.gz | wc -l
# XXX SNPs 



# --------------------
# Step 10: SAVE A COPY
# -------------------
# Save in three different locations: hard drive, cloud, etc.

# END



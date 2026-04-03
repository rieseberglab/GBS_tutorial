#####################################
#
# Part 3: SNP Calling
# March 2026
#
# Code Contributor: Yue Y & Vincent Fetterley
# Biodiversity Research Center, UBC
#
#####################################


# This code includes the following:

# Step 3.1: bwa-mem2:       fastq -> BAM
# Step 3.2: HaplotypeCaller:  BAM -> g.vcf
# Step 3.3: GenomicDBImport:         g.vcf -> database
# Step 3.4: GenotypeGVCF:                     database -> VCF
# Step 4.1: VCFtools:                                     VCF -> raw SNPs -> filtered SNPs (See script 4)




#####################################
#  3.1 bwa-mem2: fastq -> BAM
#####################################


nano run_BWA.sh
# -------------- run_BWA.sh --------------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G  # each core 8G, total 8G * 8 = 64G
#SBATCH --time=04:00:00
#SBATCH --array=1-10  # number of samples

module load StdEnv/2023
module load bwa-mem2/2.2.1
module load gcc/12.3
module load samtools/1.20

/home/yueyu/scratch/GBS_workshop/3_SNP_calling/fastq

# Define variables
REF="/home/yueyu/scratch/GBS_workshop/FASTA/praecox2.fasta"  # refernce genome FASTA

# Extract sample name from all files
i=$(ls *paired_R1.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1) 
SAMPLE=$(echo $i | cut -d "_" -f 1-2)

R1="${SAMPLE}_paired_R1.fastq.gz"
R2="${SAMPLE}_paired_R2.fastq.gz"

# OutputBAM
OUTPUT_BAM="BAM_output/${SAMPLE}.callonPRAHap2.sort.bam"

# BWA-MEM2 & samtools sort
bwa-mem2 mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" \
    "$REF" "$R1" "$R2" | \
samtools sort -@ 8 -m 7G -o "$OUTPUT_BAM"

# -- Note:
# -R: add read group (RG) information to the output BAM/SAM file. 
# It is necessary for GATK to run later on.

# -------------- run_BWA.sh (END) --------------------


# -- Index bam file --> bam.bai
cd /home/yueyu/scratch/GBS/BAM_output
module load samtools/1.20

for bam in *.sort.bam; do
    samtools index "$bam"
done




#####################################
#  3.2 HaplotypeCaller: BAM -> g.vcf
#####################################


nano run_HapCaller.sh
#----------------  run_HapCaller.sh (start)  ----------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1-10       # number of samples

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/GBS/BAM_output


# Define variables
REF="/home/yueyu/scratch/GBS_workshop/FASTA/praecox2.fasta"  # refernce genome FASTA

i=$(ls *sort.bam | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
SAMPLE=$(echo $i | cut -d "." -f 1) 

gatk --java-options "-Xmx28G" HaplotypeCaller \
	 -R ${REF} \
	 -I ${SAMPLE}.callonPRAHap2.sort.bam \
	 -O /home/yueyu/scratch/GBS/GVCF_output/${SAMPLE}.callonPRAHap2.g.vcf \
	 -ERC GVCF \
	 --max-alternate-alleles 3 \
	 --pcr-indel-model AGGRESSIVE \
	 -G StandardAnnotation -G AS_StandardAnnotation 

#----------------  run_HapCaller.sh (end)  ----------------



#####################################
#  3.3 GenomicDBImport: g.vcf -> database
#####################################


#---  Move all GVCF to the same folder
cd /home/yueyu/scratch/GBS/GVCF_output

#---  Create sample map (required to create the DBI database)
cd /home/yueyu/scratch/GBS/GVCF_output

for i in *.g.vcf; do
    name=$(echo "$i" | cut -d "." -f 1)
    echo -e "$name\t$i" >> sample_names.txt
done

head sample_names.txt


# -- Import all samples into DBI per CHROMOSOME 

#                 |------> Chr01
#  G.VCF (Genome)-|------> Chr02
#                 |------> Chr03
#                 |------> ChrXX


cd /home/yueyu/scratch/GBS/GVCF_output

for i in $(seq -w 01 17)
do
  cat <<EOL > PRA_Chr${i}_makeDB.sh
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=10G

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/GBS/GVCF_output

gatk --java-options "-Xmx48g" GenomicsDBImport \\
    --genomicsdb-workspace-path PRA_CHROM$i \\
    --batch-size 50 \\
    --sample-name-map sample_names.txt \\
    --reader-threads 3 \\
    -L PRA_chr$i
EOL

  chmod +x PRA_Chr${i}_makeDB.sh
done


# -- Make final run script
for i in {01..17};
do 
    echo "sbatch PRA_Chr${i}_makeDB.sh" >> run_DBI_byCHROM.sh
done


#-- Run in parallel 
cat run_DBI_byCHROM.sh | parallel -j 17


#-- Check status
for i in {30639408..30639424};do
    seff $i | grep State >> check_PRA_makeDB.txt
done

#-- Remove files
rm slurm-*
rm check_PRA_makeDB.txt





#####################################
#  3.4 GenotypeGVCF:  database -> VCF
#####################################

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF

for i in $(seq -w 01 17)
do
  cat <<EOL > PRA_Chr${i}_genotypeVCF.sh
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/G_VCF_ALL_merged

gatk --java-options "-Xmx90g" GenotypeGVCFs \\
     -R praecox2.fasta \\
     -V gendb://PRA_CHROM$i \\
     -O /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF/PRA_Chr$i.vcf.gz
EOL

  chmod +x PRA_Chr${i}_genotypeVCF.sh
done


#Make final run script
for i in {01..17};
do 
    echo "sbatch PRA_Chr${i}_genotypeVCF.sh" >> run_PRA_GT_byCHR.sh
done


#Run in parallel
cat run_PRA_GT_byCHR.sh | parallel -j 17





#--- Merge all CHROM into final VCF

#                 |------> Chr01
#  G.VCF (Genome)-|------> Chr02 ---> Merge, Final VCF
#                 |------> Chr03
#                 |------> ChrXX


nano run_PRA_vcf_merge.sh

#----------  vcf_merge.sh  ---------- 
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6


cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF

bcftools concat \
  PRA_Chr01.vcf.gz \
  PRA_Chr02.vcf.gz \
  PRA_Chr03.vcf.gz \
  PRA_Chr04.vcf.gz \
  PRA_Chr05.vcf.gz \
  PRA_Chr06.vcf.gz \
  PRA_Chr07.vcf.gz \
  PRA_Chr08.vcf.gz \
  PRA_Chr09.vcf.gz \
  PRA_Chr10.vcf.gz \
  PRA_Chr11.vcf.gz \
  PRA_Chr12.vcf.gz \
  PRA_Chr13.vcf.gz \
  PRA_Chr14.vcf.gz \
  PRA_Chr15.vcf.gz \
  PRA_Chr16.vcf.gz \
  PRA_Chr17.vcf.gz \
  -O z > cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF/Raw_VCF.vcf.gz

tabix -p vcf Raw_VCF.vcf.gz

#----------  vcf_merge.sh  (END) ---------- 

# -- Check number of SNPs in the raw VCF file
bcftools view -H Raw_VCF.vcf.gz | wc -l 


# END

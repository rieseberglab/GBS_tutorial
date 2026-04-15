#####################################
#
# Part 3: SNP Calling
# March 2026
#
# Code Contributor: Yue Y & Vincent F
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

cd ~/scratch/GBS_workshop
mkdir BAM_output

nano run_BWA.sh

# -------------- run_BWA.sh --------------------

#!/bin/bash
#SBATCH --account=def-prof # need to set for your PI's account
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G  # each core 8G, total 8G * 8 = 64G
#SBATCH --time=04:00:00
#SBATCH --array=1-10  # number of samples

module load StdEnv/2023
module load bwa-mem2/2.2.1
module load gcc/12.3
module load samtools/1.20

cd ~/scratch/GBS_workshop/3_SNP_calling/fastq

# Define variables
REF="/home/username/scratch/GBS_workshop/FASTA/praecox2.fasta"  # reference genome FASTA

# Extract sample name from all files
i=$(ls *paired_R1.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1) 
SAMPLE=$(echo $i | cut -d "_" -f 1-2)

R1="${SAMPLE}_paired_R1.fastq.gz"
R2="${SAMPLE}_paired_R2.fastq.gz"

# OutputBAM
OUTPUT_BAM="/home/username/scratch/GBS_workshop/BAM_output/${SAMPLE}.callonPRAHap2.sort.bam"

# BWA-MEM2 & samtools sort
#bwa-mem2 index $REF # index the reference genome

bwa-mem2 mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" \
    "$REF" "$R1" "$R2" | \
samtools sort -@ 8 -m 7G -o "$OUTPUT_BAM"

# -- Note:
# -R: add read group (RG) information to the output BAM/SAM file. 
# It is necessary for GATK to run later on.

# -------------- run_BWA.sh (END) --------------------

sbatch run_BWA.sh


# -- Index bam file --> bam.bai
cd ~/scratch/GBS_workshop/BAM_output
module load samtools/1.20

for bam in *.sort.bam; do
    samtools index "$bam"
done




#####################################
#  3.2 HaplotypeCaller: BAM -> g.vcf
#####################################

cd ~/scratch/GBS_workshop
nano run_HapCaller.sh

#----------------  run_HapCaller.sh (start)  ----------------

#!/bin/bash
#SBATCH --account=def-prof # need to set for your PI's account
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1-10       # number of samples

module load StdEnv/2020
module load gatk/4.2.4.0

mkdir ~/scratch/GBS_workshop/GVCF_output
cd ~/scratch/GBS_workshop/BAM_output


# Define variables
REF="/home/username/scratch/GBS_workshop/FASTA/praecox2.fasta"  # reference genome FASTA

i=$(ls *sort.bam | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
SAMPLE=$(echo $i | cut -d "." -f 1) 

#gatk CreateSequenceDictionary \
#    -R /home/scratch/GBS_workshop/FASTA/praecox2.fasta \
#    -O /home/username/scratch/GBS_workshop/FASTA/praecox2.dict

gatk --java-options "-Xmx28G" HaplotypeCaller \
	 -R ${REF} \
	 -I ${SAMPLE}.callonPRAHap2.sort.bam \
	 -O /home/username/scratch/GBS_workshop/GVCF_output/${SAMPLE}.callonPRAHap2.g.vcf \
	 -ERC GVCF \
	 --max-alternate-alleles 3 \
	 --pcr-indel-model AGGRESSIVE \
	 -G StandardAnnotation -G AS_StandardAnnotation 

#----------------  run_HapCaller.sh (end)  ----------------

sbatch run_HapCaller.sh





#####################################
#  3.3 GenomicDBImport: g.vcf -> database
#####################################

cd ~/scratch/GBS_workshop/GVCF_output

#---  Create sample map (required to create the DBI database)
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


cd ~/scratch/GBS_workshop/GVCF_output


#######################################
# BLOCK 1 #
#######################################


for i in $(seq -w 01 17)
do
  cat <<EOL > PRA_Chr${i}_makeDB.sh
#!/bin/bash
#SBATCH --account=def-prof # need to set for your PI's account
#SBATCH --time=10:00:00
SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=10G

module load StdEnv/2020
module load gatk/4.2.4.0

cd ~/scratch/GBS_workshop/GVCF_output

gatk --java-options "-Xmx48g" GenomicsDBImport \\
    --genomicsdb-workspace-path PRA_CHROM$i \\
    --batch-size 50 \\
    --sample-name-map sample_names.txt \\
    --reader-threads 3 \\
    -L PRA_chr$i
EOL

  chmod +x PRA_Chr${i}_makeDB.sh
done
    
#####------END BLOCK 1-----#############


# -- Make final run script

########################################
# BLOCK 2 #
########################################

for i in {01..17};
do 
    echo "sbatch PRA_Chr${i}_makeDB.sh" >> run_DBI_byCHROM.sh
done

#####------END BLOCK 2-----#############


#-- Run in parallel 
cat run_DBI_byCHROM.sh | parallel -j 17


#-- Check status
for i in {59222375..59222391};do # change according to your JOBIDs
    seff $i | grep State >> check_PRA_makeDB.txt
done
    
cat check_PRA_makeDB.txt

#-- Remove files
rm slurm-*
rm check_PRA_makeDB.txt





#####################################
#  3.4 GenotypeGVCF:  database -> VCF
#####################################

cd ~/scratch/GBS_workshop/GVCF_output

for i in $(seq -w 01 17)
do
  cat <<EOL > PRA_Chr${i}_genotypeVCF.sh
#!/bin/bash
#SBATCH --account=def-prof # need to set for your PI's account
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load gatk/4.2.4.0

cd ~/scratch/GBS_workshop/GVCF_output

gatk --java-options "-Xmx90g" GenotypeGVCFs \\
     -R ../FASTA/praecox2.fasta \\
     -V gendb://PRA_CHROM$i \\
     -O PRA_Chr$i.vcf.gz
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
#SBATCH --account=def-prof # need to set for your PI's account
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6


cd ~/scratch/GBS_workshop/GVCF_output

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
  -O z > Raw_VCF.vcf.gz

tabix -p vcf Raw_VCF.vcf.gz

#----------  vcf_merge.sh  (END) ---------- 

sbatch run_PRA_vcf_merge.sh

# -- Check number of SNPs in the raw VCF file
bcftools view -H Raw_VCF.vcf.gz | wc -l #477912


# END

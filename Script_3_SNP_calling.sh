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
#SBATCH --array=1-117  # number of samples

module load StdEnv/2023
module load bwa-mem2/2.2.1
module load gcc/12.3
module load samtools/1.20

cd /home/yueyu/scratch/GBS/fastq

# Define variables
REF="praecox2.fasta"  # refernce genome FASTA

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
#SBATCH --array=1-117       # number of samples

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/GBS/BAM_output


# Define variables
REF="praecox2.fasta"

i=$(ls *sort.bam | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
SAMPLE=$(echo $i | cut -d "." -f 1) 

gatk --java-options "-Xmx28G" HaplotypeCaller \
	 -R ${REF} \
	 -I ${SAMPLE}.callonPRAHap2.sort.bam \
	 -O GVCF_output/${SAMPLE}.callonPRAHap2.g.vcf \
	 -ERC GVCF \
	 --max-alternate-alleles 3 \
	 --pcr-indel-model AGGRESSIVE \
	 -G StandardAnnotation -G AS_StandardAnnotation 

#----------------  run_HapCaller.sh (end)  ----------------



#####################################
#  3.3 GenomicDBImport: g.vcf -> database
#####################################








#####################################
#  3.4 GenotypeGVCF:  database -> VCF
#####################################



# END

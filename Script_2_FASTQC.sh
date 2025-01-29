#####################################
#
# Part 2: fastQC
# January 2025
#
# Code Contributor: Cassandra E. & Yue Y
# Biodiversity Research Center, UBC
#
#####################################

# -- move to directory 
cd /scratch/celphin/GBS_workshop/2_fastQC

# -- load module 
module load fastqc

# -- run FASTQC
fastqc PopBARD_1.R1.fq.gz PopBARD_1.R2.fq.gz

# -- Download to own computer through Globus

# -- Visualize your fastQC results

 #END

#####################################
#
# Part 2: fastQC
# January 2025
#
# Code Contributor: Cassandra E., Yue Y & Vincent F
# Biodiversity Research Center, UBC
#
#####################################

# -- move to directory 
cd /scratch/username/GBS_workshop/2_fastQC

# -- load module 
module load fastqc

# -- run FASTQC
fastqc PopBARD_1.R1.fq.gz PopBARD_1.R2.fq.gz

# -- Download to own computer through Globus
# -- If you dont have Globus set up, you can scp to your own computer
scp username@narval.computecanada.ca:/home/username/scratch/GBS_workshop/2_fastQC/*fastqc* /path/on/own/computer

# -- Visualize your fastQC results

 #END

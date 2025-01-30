#####################################
#
# Part 0: Coding basics
# January 2025
#
# Code Contributor: Cassandra E.
# Biodiversity Research Center, UBC
#
#####################################

# Globus login: https://www.globus.org/data-transfer

# Workshop data directory: /scratch/celphin/GBS_workshop/

################################
# Terminal
# Windows users: https://mobaxterm.mobatek.net/download.html 

# Mac/Linux
ssh -X username@narval.computecanada.ca

####################################

# Exploring the server

ls # to list files in directory
cd scratch # to change directories
mkdir GBS_tutorial # to make a new directory

#----------------------------
# Tmux or screen sessions

tmux new-session -s GBS
cd scratch
cp -v /home/celphin/scratch/GBS_workshop* .
#Ctrl+B d # will exit the tmux session 

tmux attach-session -t GBS # reattach
#Ctrl+B d # will exit the tmux session again

#-------------------------------
# Modules
# https://docs.alliancecan.ca/wiki/Available_software 

module spider python
module load python
module list
python 

#------------------------------------
# Please run this before Day 2 - install R packages using tmux and module

# tmux new-session -s GBS
tmux attach-session -t GBS # reattach

module load StdEnv/2020
module load r/4.2.1

R

install.packages("tidyverse")

install.packages("ggplot2")

 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("SNPRelate")
 

#Ctrl+B d # will exit the tmux session 

#-----------------------------------------
# Quota
# https://docs.alliancecan.ca/wiki/Storage_and_file_management 

quota
               # Description            			Space       	       # of files
             	# /home (project celphin)        	2461M/50G         	54k/500k
          	# /scratch (project celphin)          	16T/20T          	86k/1000k
          	# /project (project def-prof)       	 441k/1000G        18/500k

#------------------------------------------
# Scripts

mkdir ~/scratch/test; cd ~/scratch/test
nano test_script.sh # makes and opens a script file

# can also make file with
cat << EOF > test_script.sh
#!/bin/bash
echo hello $USER
EOF

# to run
./test_script.sh

# update permissions
chmod +755 test_script.sh 

# try running again
./test_script.sh

#-------------------------------------------
# Slurm interactive allocations
# https://docs.alliancecan.ca/wiki/Running_jobs 

# tmux new-session -s GBS
tmux attach-session -t GBS # reattach
# Node: celphin@narval1

# request an interactive allocation
salloc -c1 --time 3:00:00 --mem 20000m --account def-prof
# wait for allocation to be granted
Ctrl+B d # will exit the tmux session 

sq # to view queue
seff # to see job efficiency 
exit # to cancel the allocation once you get it and are in it
scancel -u <username> # to cancel all slurm allocations or requests

#-------------------------------------------
# Slurm script example
# https://docs.alliancecan.ca/wiki/Running_jobs 

cat << EOF > slurm_script.sh
#!/bin/bash
#SBATCH --account=def-prof # need to set for your PI's account
#SBATCH --time=0-02:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G
module load StdEnv/2020
module load fastsimcoal2/2.7.0.9
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2024/runs/5pops_mertensiana_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q
EOF

# to submit the job script to the queue
sbatch slurm_script.sh 

# to check the queue
sq

#################################
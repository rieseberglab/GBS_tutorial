#####################################
#
# Part 4: ADMIXTURE 
# January 2025
#
# Code Contributor: Cassandra E. & Yue Y.
# Biodiversity Research Center, UBC
#
#####################################

# This code includes the following:
# 4.1: make BED files from filtered VCF files
# 4.2: run admixture for different k values with loop
# 4.3: bootstrap admixture results
# 4.4: Plot ADMIXTURE result in R



#####################################
#  4.1 VCF -> BED
#####################################
# make Plink file to be able to run Admixture

cd /scratch/celphin/GBS_workshop/4_ADMIXTURE

# Copy the filtered VCF file from STEP 3 to the current directory
cp /scratch/celphin/GBS_workshop/3_SNP_filtering/Cassiope_noMER_r10i.recode.vcf .

# convert the vcf file to Plink
vcftools --vcf Cassiope_noMER_r10i.recode.vcf --plink --out Cassiope_noMER_r10i

#MAKE A BED FILE 
plink --file Cassiope_noMER_r10i --allow-no-sex --allow-extra-chr 0 --make-bed --out Cassiope_noMER_r10i

#Save sample names
module load bcftools
bcftools query -l Cassiope_noMER_r10i.recode.vcf > sample_39_names.txt


#####################################
#  4.2 run ADMIXTURE
#####################################
# see https://dalexander.github.io/admixture/admixture-manual.pdf for more details

tmux new-session -s Cassiope
tmux attach-session -t Cassiope
#salloc -c48 --time 03:00:00 --mem 30G --account def-rieseber

cd /scratch/celphin/GBS_workshop/4_ADMIXTURE

module load StdEnv/2020
module load admixture/1.3.0

mkdir Take1
mkdir bootstrap

cp Cassiope_noMER_r10i* Take1
cp Cassiope_noMER_r10i* bootstrap/

cd Take1
for K in 1 2 3 4 5 6 7 8 9; \
do admixture --cv=10 -s time -j48 -C 0.0000000001  Cassiope_noMER_r10i.bed $K | tee log${K}.out; done


#to get the CV errors and see which K value is the best model
grep -h CV log*out

# Take 1
CV error (K=1): 0.54249
*CV error (K=2): 0.33333
CV error (K=3): 0.34385
CV error (K=4): 0.37351
CV error (K=5): 0.37667
CV error (K=6): 0.42367



#####################################
#  4.3 run ADMIXTURE (PERMUTATION 1000) - DO NOT RUN FOR THE WORKSHOP
#####################################

cd /scratch/celphin/GBS_workshop/4_ADMIXTURE/bootstrap

module load StdEnv/2020
module load admixture/1.3.0

#make one last run of K=2-6 with 1000 bootstraps to estimate Standard Errors of the estimates

for K in 2 3 4 5 6 ; \
do admixture -B1000 --cv=10 -s time -j48 -C 0.0000000001 Cassiope_noMER_r10i.bed $K | tee log${K}.out; done

# Yue Note: this is probably going to take a long time, are we running this in the workshop to check standarded errors




#####################################
#  4.4: Plot ADMIXTURE in R
#####################################

cd /scratch/celphin/GBS_workshop/4_ADMIXTURE

cp sample_40_names.txt /Take1

module load StdEnv/2020
module load r/4.2.1

R

# install.packages("tidyverse")
# install.packages("ggplot2")
# pkg install takes a while (~15min)

library(tidyverse)
library(ggplot2)


# set the workimg directory
getwd()
setwd("/scratch/celphin/GBS_workshop/4_ADMIXTURE/Take1")

samplelist <- read_tsv("sample_39_names.txt",col_names = "sample")

all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())


for (k in 1:6){
  data <- read_delim(paste0("Cassiope_noMER_r10i.",k,".Q"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")

  data$sample <- samplelist$sample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)

  }

all_data


#Plot

plot <- all_data %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=seq(1:5)) +
  facet_wrap(~k,ncol=1)


#Change width and length!! (according to need)
png(filename = "Cassiope.admixture.k1to6.png", width = 1400, height = 1100, res = 200) 
print(plot)
dev.off()


# Using Globus, move to local computer, check final plot


# ---------- END (Yue 2025 Jan 28th)


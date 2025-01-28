
#####################################
#
# Part 4: PCA & Plotting in R
# January 2025
#
# Code Contributor: Cassandra E.
# Biodiversity Research Center, UBC
#
#####################################


# This code includes the following:
# 4.1: Load R pkgs
# 4.2: Load data for plotting
# 4.3: Run & Plot PCA





# ------------- TO - DO --------------- 

#####################################
#  4.1 Load R pkgs
#####################################

cd /scratch/celphin/GBS_workshop/4_PCA

# open R in terminal

module load StdEnv/2020
module load r/4.0.0

R

# -- install packages
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("SNPRelate")

# -- load pkgs
library(tidyverse)
library(ggplot2)
library(SNPRelate)


#####################################
#  4.2 Load data for plotting
#####################################


setwd("~/GitHub/Population_genomics_Cassiope")

#Load the gds file - for PCA
genofile <- snpgdsOpen("./Figures_data/Cassiope_noMER_r10i.recode.gds")

# Sample names 40 in filtered VCF
samples_list <- read.table("./Figures_data/Cassiope_noMER_r10i.imiss", header = TRUE)



#####################################
#  4.3 Run & Plot PCA
#####################################



#################################
# PCA
# https://owensgl.github.io/biol525D/Topic_8-9/pca.html

# create GDS file on server
# snpgdsVCF2GDS("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.vcf.gz",
# "/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.gds",
# method="biallelic.only")

#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

#Make a list of sites we're keeping.
snpset.id <- unlist(snpset_pruned)

#Run the PCA
pca0 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, autosome.only = F)
pca5 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, maf=0.05,  autosome.only = F)

###################################
pca <- pca5

#Here's the percent variance explained for each eigenvector
pc.percent <- pca$varprop*100
round(pc.percent, 2)
# MAF 1%
# 27.00  5.62  2.63  1.63  1.50  1.26  1.02  0.93  0.93  0.85  0.82  0.77  0.73  0.73  0.68  0.67
# MAF 5%
# 10.08  5.25  3.57  3.28  1.86  1.57  1.37  1.26  1.20  1.19  1.10  0.92  0.82  0.79  0.78  0.76

#Make a dataframe of your PCA results
PCA_tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  PC3 = pca$eigenvect[,3],    # the first eigenvector
                  PC4 = pca$eigenvect[,4],    # the second eigenvector
                  PC5 = pca$eigenvect[,5],    # the first eigenvector
                  PC6 = pca$eigenvect[,6],    # the second eigenvector
                  stringsAsFactors = FALSE)

str_split(as.character(PCA_tab$sample), "_")

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(PCA_tab$sample), "_"))))))
PCA_tab_data <- cbind(PCA_tab, q[,1])

colnames(PCA_tab_data) <- c("ID_code", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6","Pop")


# PCA DONE








#####################################

#Plot a PCA image
# MAF 1%
# PCA = 27.00  5.62  2.63  1.63  1.50  1.26  1.02  0.93  0.93  0.85  0.82  0.77  0.73  0.73  0.68  0.67

colnames(All_samples_data)[2:6] <- c("Russia", "Europe", "Subspp", "Greenland", "Alaska")

png("./Figures_data/Plots/PCA_PC1_PC2_maf1.png", width = 3000, height = 2700)
All_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = Max_Admix_Group), size=15)  +
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC2 (5.62%)", x = "PC1 (27.00%)")
dev.off()

png("./Figures_data/Plots/PCA_PC2_PC3_maf1.png", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC2,y=PC3)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC3 (2.63%)", x = "PC2 (5.62%)")
dev.off()

png("./Figures_data/Plots/PCA_PC3_PC4_maf1.png", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC3,y=PC4)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC4 (1.63%)", x = "PC3 (2.63%)")
dev.off()

png("./Figures_data/Plots/PCA_PC5_PC6_maf1.png", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC5,y=PC6)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC6 (1.26%)", x = "PC5 (1.50%)")
dev.off()

colnames(All_samples_data)[2:6] <- c("V1", "V2", "V3", "V4", "V5")






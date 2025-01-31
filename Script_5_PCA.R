
#####################################
#
# Part 5: PCA & Plotting in R
# January 2025
#
# Code Contributor: Cassandra E. & Yue Y.
# Biodiversity Research Center, UBC
#
#####################################


# This code includes the following:
# 5.1: Load R pkgs
# 5.2: Load data for plotting
# 5.3: Run PCA
# 5.4: Plot PCA



#####################################
#  5.1 Load R pkgs
#####################################

cd /scratch/celphin/GBS_workshop/5_PCA

# -- open R in terminal

module load StdEnv/2020
module load r/4.2.1

R

# -- install packages
# install.packages("tidyverse")
# install.packages("ggplot2")
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("SNPRelate")

# -- load pkgs
library(tidyverse)
library(ggplot2)
library(SNPRelate)


#####################################
#  5.2 Load data for plotting
#####################################

# -- Change this to your folder pathway 
setwd("/scratch/celphin/GBS_workshop/5_PCA")

# -- Load the gds file - for PCA


# -- Sample names 40 in filtered VCF
samples_list <- read.table("/scratch/celphin/GBS_workshop/4_ADMIXTURE/sample_39_names.txt", header = TRUE)


# -- VCF to GDS (IMPORTANT STEP)
#    SNPRelate works with a compressed version of a genotype file called a “gds”
snpgdsVCF2GDS("/scratch/celphin/GBS_workshop/4_ADMIXTURE/Cassiope_noMER_r10i.recode.vcf",
              "Cassiope_noMER_r10i.recode.gds",
              method="biallelic.only")


# -- Then load GDS file
genofile <- snpgdsOpen("/scratch/celphin/GBS_workshop/5_PCA/Cassiope_noMER_r10i.recode.gds")


# -- Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

# -- Make a list of sites we're keeping.
snpset.id <- unlist(snpset_pruned)



#####################################
#  5.3 Run PCA
#####################################


# -- Run the PCA
pca <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, autosome.only = F)


# -- Here's the percent variance explained for each eigenvector
pc.percent <- pca$varprop*100
round(pc.percent, 2)
# 30.58  6.02  3.85  2.92  2.75  2.71  2.58  2.54  2.47  2.29  2.26  2.24 2.21  2.12  2.08  2.04


# -- Make a dataframe of your PCA results
PCA_tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the 1st eigenvector
                  PC2 = pca$eigenvect[,2],    # the 2nd eigenvector
                  PC3 = pca$eigenvect[,3],    # the 3rd eigenvector
                  PC4 = pca$eigenvect[,4],    # the 4th eigenvector
                  PC5 = pca$eigenvect[,5],    # the 5th eigenvector
                  PC6 = pca$eigenvect[,6],    # the 6th eigenvector
                  stringsAsFactors = FALSE)

head(PCA_tab)
dim(PCA_tab)

population_ID <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(PCA_tab$sample), "_"))))))
PCA_tab_final <- cbind(PCA_tab, population_ID[,1])
colnames(PCA_tab_final) <- c("ID_code", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6","Pop")
PCA_tab_final$Pop <- gsub("[0-9]", "", PCA_tab_final$Pop)

head(PCA_tab_final)
dim(PCA_tab_final)

unique(PCA_tab_final$Pop)




#####################################
#  5.4 Plot PCA
#####################################


# -- Format to plot in command line ---
# png("filename.png",width = XX, height = XX)
# plot()
# dev.off()
# -- Format to plot in command line (END) ---


# PC1    PC2   PC3   PC4   PC5   PC6 
# 30.58  6.02  3.85  2.92  2.75  2.71 


# -------- Plot: PC 1 V.S. PC2 ---------
png("PCA_PC1_PC2.png", width = 1000, height = 800)
PCA_tab_final %>%
  ggplot(.,aes(x = PC1,y = PC2, color = Pop)) +
  geom_point(size = 10)  +
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 30),
        axis.text = element_text(size=30),
        axis.title=element_text(size=30)) +
  labs(y= "PC2 (6.02%)", x = "PC1 (30.58%)")
dev.off()



# -------- Plot: PC 2 V.S. PC3 ---------
png("PCA_PC2_PC3.png", width = 1000, height = 800)

PCA_tab_final %>%
  ggplot(.,aes(x=PC2,y=PC3, color = Pop)) + 
  geom_point(size=10)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 30),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30)) + 
  labs(y= "PC3 (3.85%)", x = "PC2 (6.02%)")
dev.off()




# -------- Plot: PC 3 V.S. PC4 ---------
png("PCA_PC3_PC4.png", width = 1000, height = 800)

PCA_tab_final %>%
  ggplot(.,aes(x=PC3,y=PC4, color = Pop)) + 
  geom_point(size=10)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 30),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30)) + 
  labs(y= "PC4 (2.92%)", x = "PC3 (3.85%)")
dev.off()



# -------- Plot: PC 5 V.S. PC6 ---------
png("PCA_PC5_PC6.png", width = 1000, height = 800)

PCA_tab_final %>%
  ggplot(.,aes(x=PC5,y=PC6, color = Pop)) + 
  geom_point(size=10)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 30),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30)) + 
  labs(y= "PC6 (2.71%)", x = "PC5 (2.75%)")
dev.off()






# ---------- END (Yue 2025 Jan 28th)

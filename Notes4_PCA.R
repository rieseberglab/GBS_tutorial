
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
# 4.3: Plot ADMIXTURE results


# 2025 Jan 27th 4pm: Yue Notes
# Do we need all the packages?? it will take a while for people to install all the packages


module load StdEnv/2020
module load r/4.2.1

R

library(tidyverse)
library(ggplot2)

































# ------------- TO - Review with Cssandra her code --------------- 

#####################################
#  4.1 Load R pkgs
#####################################

cd /scratch/celphin/GBS_workshop/4_PCA

# open R in terminal

module load StdEnv/2020
module load r/4.0.0

R

# -- install packages
 install.packages("devtools") 


 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("SNPRelate")
 
 install.packages(c("tidyverse","dplyr","vctrs", "tidyr","reshape","geosphere","ade4","ape","scatterpie","data.table","maps","readr","stringr","ggplot2"))

# -- load pkgs
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape)
library(geosphere)
library(ade4)
library(ape)
library(scatterpie)
library(data.table)
library(maps)
library(readr) 
library(stringr)
library(SNPRelate)
library(ggplot2)



# Not


#####################################
#  4.2 Load data for plotting
#####################################


setwd("~/GitHub/Population_genomics_Cassiope")

#Load the gds file - for PCA
genofile <- snpgdsOpen("./Figures_data/Cassiope_noMER_r10i.recode.gds")

# load list of samples in filtered vcf
samples_list <- read.table("./Figures_data/Cassiope_noMER_r10i.imiss", header = TRUE)

# load in Admixture data - code setup for 5 populations
Admix_tbl=read.table("./Figures_data/Admixture/Take1/Cassiope_noMER_r10i.5.Q")

# load detailed information about all 371 individuals
sample_information  <- read.csv("./Figures_data/All_samples_lat_long_leaf_icetimes.csv", header = TRUE)

# load in population specific information about sites
population_information  <- read.csv("./Figures_data/Pop_LatLong_Icetimes_Climate.csv", header = TRUE)

# PopStats
Het_data <- read.table("./Figures_data/PopStats/Het_data_by_pop.txt", header = TRUE)
Pi_wind_data <- read.table("./Figures_data/PopStats/summed_wind_Pi.txt", header=TRUE, sep=" ")
Pi_site_data <- read.table("./Figures_data/PopStats/summed_site_Pi.txt", header=TRUE, sep=" ")
Tajima_data <- read.table("./Figures_data/PopStats/population_TajimaD.txt", header=TRUE, sep=" ")















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













#####################################
# Join all the data

# Sample specific information

# samples_list 
# Admix_tbl
# sample_information 
# Het_data
# PCA_tab_data

# Join Samples information
samples_list$ID_code <- as.factor(samples_list$INDV)
ID_code = samples_list$ID_code
sample_information$Pop[which(sample_information$Pop=="25")] <- "025"

Admix_tbl1 <- cbind(ID_code, Admix_tbl)
All_samples_data <- left_join(Admix_tbl1, sample_information, by="ID_code")
All_samples_data <- left_join(All_samples_data, Het_data, by="ID_code")
All_samples_data <- left_join(All_samples_data, samples_list, by="ID_code")
All_samples_data <- left_join(All_samples_data, PCA_tab_data, by="ID_code")

#------------------------------
# join all Population information

# population_information  
# Pi_wind_data 
# Pi_site_data 
# Tajima_data 

colnames(Pi_site_data) <- c("Pop" ,   "SumSitePi" , "Site_SD_Pi",  "Site_var_Pi")
colnames(Pi_wind_data) <- c("Pop" ,   "SumWindPi",  "Wind_SD_Pi",  "Wind_var_Pi")
colnames(Tajima_data ) <- c("Pop" ,   "TajimaAvg",  "TajimaSD",  "TajimaVar")
population_information$Pop[which(population_information$Pop=="25")] <- "025"

All_pop_data <- left_join(population_information, Tajima_data, by="Pop")
All_pop_data <- left_join(All_pop_data, Pi_site_data, by="Pop")
All_pop_data <- left_join(All_pop_data, Pi_wind_data, by="Pop")

#------------------------------
# average the sample information by population
Pop_avg <- All_samples_data  %>%
  group_by(Pop.x) %>%
  dplyr::summarise(V1 = mean(V1, na.rm=TRUE), 
                   V2 = mean(V2, na.rm=TRUE), 
                   V3 = mean(V3, na.rm=TRUE), 
                   V4 = mean(V4, na.rm=TRUE), 
                   V5 = mean(V5, na.rm=TRUE), 
                   # V6 = mean(V6, na.rm=TRUE), 
                   # V7 = mean(V7, na.rm=TRUE), 
                   # V8 = mean(V8, na.rm=TRUE), 
                   # V9 = mean(V9, na.rm=TRUE), 
                   # V10 = mean(V10, na.rm=TRUE),
                   Lat = mean(Lat, na.rm=TRUE),
                   Long = mean(Long, na.rm=TRUE),
                   Lat_shift = mean(Lat_shift, na.rm=TRUE),
                   Long_shift = mean(Long_shift, na.rm=TRUE),
                   Avg.Temp = mean(Avg.Temp, na.rm=TRUE),
                   Leaf_weight = mean(Leaf_weight, na.rm=TRUE),
                   Ice_time = mean(Ice_time, na.rm=TRUE),
                   O.HOM. = mean(O.HOM., na.rm=TRUE),
                   E.HOM. = mean(E.HOM., na.rm=TRUE),
                   N_SITES = mean(N_SITES, na.rm=TRUE),
                   FIS = mean(F, na.rm=TRUE),
                   N_MISS = mean(N_MISS, na.rm=TRUE),
                   F_MISS = mean(F_MISS, na.rm=TRUE),
                   PC1 = mean(PC1, na.rm=TRUE),
                   PC2 = mean(PC2, na.rm=TRUE),
                   PC3 = mean(PC3, na.rm=TRUE),
                   PC4 = mean(PC4, na.rm=TRUE),
                   PC5 = mean(PC5, na.rm=TRUE),
                   PC6 = mean(PC6, na.rm=TRUE)
                   )

Pop_avg$Pop <- Pop_avg$Pop.x
All_pop_data <- left_join(All_pop_data, Pop_avg, by="Pop")


############################################
# calculate the max Admix group for each location
Admix_groups <- gather(All_pop_data, Group, Amount, V1:V5, factor_key=TRUE)

maxAdmixgroup <- Admix_groups  %>%
  group_by(Pop) %>%
  dplyr::summarise(Amount = max(Amount, na.rm=TRUE))

Final_AdmixGroups <- left_join(maxAdmixgroup, Admix_groups, by="Amount")

All_pop_data$Max_Admix_Group <- Final_AdmixGroups$Group

Admix_K_Region <- All_pop_data  %>%
  group_by(Group) %>%
  dplyr::summarise(Max_Admix_Group = list(Max_Admix_Group))

# 5 groups
# V5 - Alaska/NWT
# V2 - Europe
# V4 - Greenland
# V1 - Russia
# V3 - Saximontana

Max_Admix_Group <- c("V1", "V2", "V3", "V4", "V5")
Admix_Location <- c("Russia", "Europe", "Saximontana", "Greenland", "Alaska")
Group_location <- as.data.frame(cbind(Max_Admix_Group, Admix_Location))

All_pop_data <- left_join(All_pop_data, Group_location, by="Max_Admix_Group")

#5 groups
Groups_summary_admix <- All_pop_data %>%
  group_by(Admix_Location) %>%
  dplyr::summarise(Locations = list(Pop), 
                   Lat = mean(Lat.x, na.rm=TRUE),
                   Long = mean(Long.x, na.rm=TRUE),
                   V1 = mean(V1, na.rm=TRUE), 
                   V2 = mean(V2, na.rm=TRUE), 
                   V3 = mean(V3, na.rm=TRUE), 
                   V4 = mean(V4, na.rm=TRUE), 
                   V5 = mean(V5, na.rm=TRUE))


###################################
# Set up colour scheme 

map_colours_5g <- c("red", "yellow",  "purple", "orange", "deepskyblue")

# Alaska - "deepskyblue"
# Greenland - "lightorange"  
# Saximontana - "purple"
# Europe - "yellow"
# Russia - "red"

map_col <- cbind(Max_Admix_Group, map_colours_5g)
map_col <- as.data.frame(map_col)

All_pop_data <- left_join(All_pop_data, map_col, by="Max_Admix_Group")

# join map_colours with individual sample data
Map_col_pop <- as.data.frame(cbind(All_pop_data$Pop, All_pop_data$map_colours_5g, All_pop_data$Max_Admix_Group))

colnames(Map_col_pop) <- c("Pop", "map_colours_5g", "Max_Admix_Group")
All_samples_data$Pop <- All_samples_data$Pop.x

All_samples_data <- left_join(All_samples_data, Map_col_pop, by="Pop")

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

#########################################
# Map of samples
# https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html

All_pop_data  <- transform_coord(All_pop_data, lon = "Long.x", lat = "Lat.x", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

All_pop_data$Pop <- as.factor(All_pop_data$Pop)

png("./Figures_data/Plots/Cassiope_samples_map_colours.png", width = 2000, height = 1400)
basemap("panarctic", limits=50) + 
  geom_point(aes(x = lon.utm, y = lat.utm), data = All_pop_data, col=All_pop_data$map_colours_5g, size=15) 
dev.off()

############################################

#Make Admixture barplot

#Excel: =RIGHT(A2,LEN(A2) - SEARCH("_", A2, SEARCH("_", A2) + 1))

mergedAdmixTable <- All_samples_data[,c("Pop.x", "ID_code","V1", "V2","V3", "V4" ,"V5", "Lat", "Long")]

ordered1 = mergedAdmixTable[order(mergedAdmixTable$Long),]
rownames(mergedAdmixTable) <- mergedAdmixTable$ID_code

row_order <- c( "PopGEN_10", "PopGEN_1",  "PopGEN_2",  "PopGEN_3",  "PopGEN_5",  "PopGEN_6", 
                "PopGEN_7",  "PopGEN_8",  "PopGEN_9",
                "PopKL1_7", "PopKL1_9",  "PopKL4_1",  "PopKL4_4",  "PopKL3_1",  "PopKL3_7",  "PopKL5_1",  "PopKL5_7",  "PopKL2_1", 
                "PopKL2_2",  "PopPC5_1",  "PopPC5_2",  "PopPC4_1",  "PopPC4_2",  "PopPC3_1",  "PopPC3_2",  "PopPC2_1", 
                "PopPC2_2",  "PopPC1_2",  
                "PopATQ_17", "PopATQ_18", "PopATQ_19", "PopATQ_21", "PopATQ_22", "PopATQ_23", "PopATQ_25", "PopATQ_30",    
                "PopATQ_36", "PopATQ_37", "PopBARD_10",    "PopBARD_11",    "PopBARD_14",    "PopBARD_15",    "PopBARD_16","PopBARD_1",    
                "PopBARD_2", "PopBARD_37",    "PopBARD_38",    "PopBARD_9", "PopMAT_11", "PopMAT_13", "PopMAT_18", "PopMAT_1", 
                "PopMAT_2",  "PopMAT_3",  "PopMAT_5",  "PopMAT_7",  "PopMAT_8",  "PopMNT_12", "PopMNT_17", "PopMNT_18",    
                "PopMNT_21", "PopMNT_22", "PopMNT_23", "PopMNT_4",  "PopMNT_5",  "PopMNT_7",  "PopMNT_9",  "PopIMN_11",    
                "PopIMN_13", "PopIMN_15", "PopIMN_17", "PopIMN_1",  "PopIMN_20", "PopIMN_21", "PopIMN_2",  "PopIMN_5", 
                "PopIMN_7",  "PopSAG_13", "PopSAG_14", "PopSAG_16", "PopSAG_21", "PopSAG_23", "PopSAG_24", "PopSAG_32",    
                "PopSAG_38", "PopSAG_5",  "PopSAG_9",  "PopDEN_15", "PopDEN_1",  "PopDEN_2",  "PopDEN_2u", "PopDEN_3", 
                "PopDEN_3u", "PopDEN_4",  "PopDEN_4u", "PopDEN_55", "PopDEN_5",  "PopMIL_11", "PopMIL_14", "PopMIL_17",    
                "PopMIL_1",  "PopMIL_20", "PopMIL_25", "PopMIL_2",  "PopMIL_37", "PopMIL_5",  "PopMIL_6",  "PopQHI_38",    
                "PopQHI_37", "PopQHI_34", "PopQHI_33", "PopQHI_32", "PopQHI_31", "PopQHI_24", "PopQHI_22", 
                "PopPEA_11", "PopPEA_18", "PopPEA_25", "PopPEA_26", "PopPEA_27",    
                "PopPEA_28", "PopPEA_30", "PopPEA_34", "PopPEA_36", "PopKUQ_14", "PopKUQ_20", "PopKUQ_26", "PopKUQ_37",    
                "PopKUQ_4",  "PopKUQ_9",  "PopYAM_10", "PopYAM_1",  "PopYAM_3",  "PopYAM_4",  "PopYAM_5",  "PopYAM_6", 
                "PopYAM_7",  "PopYAM_8",  "PopYAM_9",  "PopAXE_11", "PopAXE_13", "PopAXE_14", "PopAXE_21", "PopAXE_26",    
                "PopAXE_27", "PopAXE_2",  "PopAXE_34", "PopAXE_6",     
                "PopCR_10",  "PopCR_12",  "PopCR_13",  "PopCR_17",  "PopCR_1",   "PopCR_2",   "PopCR_4",   "PopCR_5",  
                "PopCR_7",   "PopCR_9",   "PopIq_10", "PopIq_13",  "PopIq_15",  "PopIq_16",  "PopIq_17",  "PopIq_20", 
                "PopIq_2",   "PopIq_3",   "PopIq_9","PopKik_11", "PopKik_13", "PopKik_16", "PopKik_17", "PopKik_18", "PopKik_1", 
                "PopKik_20", "PopKik_2",  "PopKik_3",  
                "PopDLG_10", "PopDLG_11", "PopDLG_12", "PopDLG_13", "PopDLG_14",    
                "PopDLG_16", "PopDLG_17", "PopDLG_7",  "PopDLG_8",  "PopDLG_9",  "PopDQG_10", "PopDQG_11", "PopDQG_12",    
                "PopDQG_13", "PopDQG_14", "PopDQG_1",  "PopDQG_2",  "PopDQG_3",  "PopDQG_4",  "PopDQG_6",  "Pop025_10",    
                "Pop025_11", "Pop025_12", "Pop025_15", "Pop025_20", "Pop025_24", "Pop025_5",  "Pop025_6",  "Pop025_7", 
                "Pop025_9",  "PopIG_10",  "PopIG_12",  "PopIG_14",  "PopIG_2",   "PopIG_3",   "PopIG_4",   "PopIG_5",  
                "PopIG_6",   "PopIG_9",   "PopEUR_1",  "PopEUR_22", "PopEUR_34", "PopEUR_5", 
                "PopFOS_2",  "PopFOS_3",  "PopFOS_7",  "PopFOS_8",  "PopHAZ_18", "PopHAZ_20", "PopHAZ_22",    
                "PopHAZ_24", "PopHAZ_26", "PopHAZ_29", "PopHAZ_33", "PopHAZ_34", "PopHAZ_6",  "PopHAZ_9",
                "PopGF_1",   "PopGF_4",   "PopSVO_38", "PopBY1_14",    
                "PopBY1_1",  "PopBY1_5",  "PopBY1_7",  "PopBY2_12", "PopBY2_3",  "PopBY2_9",  "PopBY3_14", "PopBY3_18",    
                "PopBY3_19", "PopAlexNew_49", "PopAlexNew_51", "PopAlexNew_46", "PopAlexNew_47", "PopAlexNew_36", 
                "PopAlexNew_28", "PopAlexNew_24", "PopAlexNew_6",  "PopAlexNew_21", "PopAlexNew_68", "PopAlexNew_73", 
                "PopAlexNew_33", "PopAlexNew_32", "PopAlexNew_42","PopAlexOld_37", "PopAlexOld_24","PopAlexOld_89", 
                "PopAlexOld_39", "PopAlexOld_30", "PopAlexOld_33",  
                "PopZAC_1",  "PopZAC_20", "PopZAC_21", "PopZAC_22", "PopZAC_23", "PopZAC_24",    
                "PopZAC_25", "PopZAC_26", "PopZAC_28", "PopZAC_7",  "PopLON_11", "PopLON_14", "PopLON_15", "PopLON_17",    
                "PopLON_19", "PopLON_23", "PopLON_28", "PopLON_36", "PopLON_4",  "PopLON_9",  "PopPET_1",  "PopPET_2", 
                "PopPET_3",  "PopPET_4",  "PopPET_5",  "PopLAJ_13", "PopLAJ_16", "PopLAJ_17", "PopLAJ_19", "PopLAJ_20",    
                "PopLAJ_21", "PopLAJ_23", "PopLAJ_25", "PopLAJ_26", "PopLAJ_5",  "PopSW_18",  "PopSW_22",  "PopSW_29", 
                "PopSW_2",   "PopSW_31",  "PopSW_34",  "PopSW_36",  "PopSW_3",   "PopSW_45",  "PopSW_5",   "PopYED_1a",    
                "PopYED_1",  "PopYED_2a", "PopYED_2",  "PopYED_3a", "PopYED_3",  "PopYED_4a", "PopYED_4",  "PopYED_5a",    
                "PopYED_5",  "PopSAM_1a", "PopSAM_1",  "PopSAM_2a", "PopSAM_2",  "PopSAM_3a", "PopSAM_3",  "PopSAM_4", 
                "PopSAM_5a", "PopSAM_5" )

ordered_pop <- mergedAdmixTable[match(row_order, mergedAdmixTable$ID_code), ]

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

#----------------------------------
#5 Groups

ordered_pop$Pop1 <-  as.factor(ordered_pop$Pop.x)

png("./Figures_data/Plots/Admix_5_bar.png", width = 6000, height = 3000)
# bottom, left, top, and right
par(mar=c(30,15,4,2))
barplot(t(as.matrix(ordered_pop[,c(3:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(ordered_pop$Pop1), las=2, cex.names=8, cex.axis=8)
dev.off()

####################################
# Plot Ancient and present day ancestry proportions

#Alex Old
Ancient_New_Alex <- ordered_pop[which(ordered_pop$Pop1=="AlexNew"|ordered_pop$Pop1=="AlexOld"),]

Ancient_New_Alex$ID_code <- as.factor(Ancient_New_Alex$ID_code)

orderedAncient = Ancient_New_Alex[order(Ancient_New_Alex$ID_code),]

orderedAncient$V99 <- orderedAncient$V10

png("./Figures_data/Plots/Admix_Ancient_bar5.png", width = 1000, height = 500)
barplot(t(as.matrix(orderedAncient[,c(3:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(orderedAncient$ID_code), las=2, cex.names=1.4)
dev.off()

#---------------------
# test if significant difference in proportions
# https://www.r-tutor.com/elementary-statistics/non-parametric-methods/mann-whitney-wilcoxon-test

# for Ellesmere (yellow)
wilcox.test(V2 ~ Pop.x, data=Ancient_New_Alex) 
#W = 7, p-value = 0.00227
#alternative hypothesis: true location shift is not equal to 0

# for Alaska (blue)
wilcox.test(V5 ~ Pop.x, data=Ancient_New_Alex) 
#W = 75, p-value = 0.004644
#alternative hypothesis: true location shift is not equal to 0

#-------------------
# barplot of difference for all V1-V5

gath_ancient <- gather(orderedAncient,"popGroup", "prob", V1:V5)

png("./Figures_data/Plots/Ancient_Admix_Amounts.png", width = 500, height = 500)
ggplot(data=gath_ancient, aes(x=Pop.x, y=prob))+
  geom_boxplot(aes(x=Pop.x, y=prob, fill="blue"))+
  facet_wrap(~popGroup, scales = "free")+ 
  theme_classic()+
  scale_fill_manual(values="blue")
  #labs(fill="popGroup")
dev.off()

# plot V2

png("./Figures_data/Plots/Ancient_V2_Amount.png", width = 400, height = 300)
ggplot(data=orderedAncient, aes(x=Pop.x, y=V2))+
  geom_boxplot(aes(x=Pop.x, y=V2, fill="yellow"))+
  theme_classic()+
  labs(y= "Europe/Ellesmere ancestry", x = "")+
  theme(axis.text = element_text(size = 20), text = element_text(size = 20)) +
  scale_fill_manual(values="yellow")
#labs(fill="popGroup")
dev.off()


# plot V5

png("./Figures_data/Plots/Ancient_V5_Amount.png", width = 400, height = 300)
ggplot(data=orderedAncient, aes(x=Pop.x, y=V5), cex=6)+
  geom_boxplot(aes(x=Pop.x, y=V5, fill="deepskyblue"))+
  theme_classic()+
  labs(y= "Alaska ancestry", x = "")+
  theme(axis.text = element_text(size = 20), text = element_text(size = 20)) +
  scale_fill_manual(values="deepskyblue")
#labs(fill="popGroup")
dev.off()


#------------------------
#Ellesmere
Ancient_New_Alex <- ordered_pop[which(ordered_pop$Pop1=="AlexNew"|ordered_pop$Pop1=="AlexOld"|ordered_pop$Pop1=="EUR"|ordered_pop$Pop1=="AXE"|ordered_pop$Pop1=="PopFOS"|ordered_pop$Pop1=="GF"|ordered_pop$Pop1=="HAZ"|ordered_pop$Pop1=="SVN"|ordered_pop$Pop1=="SVO"),]
Ancient_New_Alex$ID_code <- as.factor(Ancient_New_Alex$ID_code)
orderedAncient = Ancient_New_Alex[order(Ancient_New_Alex$ID_code),]

png("./Figures_data/Plots/Admix_Ellesmere_bar5.png", width = 1000, height = 707)
barplot(t(as.matrix(orderedAncient[,c(3:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(orderedAncient$ID_code), las=2)
dev.off()


#####################################
# Admixture map - shifted to see all groups

#library(devtools)
#devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
#library(PlotSvalbard)
# https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html
#library(scatterpie)

# shifted coordinates
All_samples_data <- transform_coord(All_samples_data, lon = "Long_shift", lat = "Lat_shift", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")
All_samples_data$lon_shift.utm <- All_samples_data$lon.utm 
All_samples_data$lat_shift.utm <- All_samples_data$lat.utm 
All_samples_data <- subset(All_samples_data, select = -c(lon.utm, lat.utm))

# regular coordinates
All_samples_data <- transform_coord(All_samples_data, lon = "Long", lat = "Lat", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

# make dataset of just what needs to be plotted
shift_Admix_lat_long <- subset(All_samples_data, select = c(ID_code, lon.utm, lat.utm, lon_shift.utm, lat_shift.utm, Pop, V1, V2, V3, V4, V5))

Pop <- c("V1", "V2", "V3", "V4", "V5")

png("./Figures_data/Plots/Admix_map5_take1_shiftedlat_long.png", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon_shift.utm, y = lat_shift.utm, group = ID_code, r=200000), data = shift_Admix_lat_long, cols = Pop, size = 0.9) +
  geom_point(aes(x = lon_shift.utm, y = lat_shift.utm), data = shift_Admix_lat_long, col="white", size=17) +
  geom_point(aes(x = lon.utm, y = lat.utm), data = shift_Admix_lat_long, col="green", size=6) +
  scale_fill_manual(values=map_colours_5g)+
  geom_text(data = shift_Admix_lat_long, aes(x = lon_shift.utm, y = lat_shift.utm, label = Pop), color = "black", fontface = 2, size = 25.4/72.27*20)
dev.off()

#---------------------------------
# Averaged admixture map - not plotting each individual

Pop <- c("V1", "V2", "V3", "V4", "V5")

# map showing population structure makeup
png("./Figures_data/Plots/Admix_map5g_take1_avg.png", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon.utm, y = lat.utm, group = Pop, r=200000), data = All_pop_data, cols = Pop, size = 0.5) +
  scale_fill_manual(values=map_colours_5g)
dev.off()

#----------------------------
All_pop_data <- subset(All_pop_data, select = -c(lon.utm, lat.utm))

All_pop_data <- transform_coord(All_pop_data, lon = "Long_shift", lat = "Lat_shift", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")
All_pop_data$lon_shift.utm <- All_pop_data$lon.utm 
All_pop_data$lat_shift.utm <- All_pop_data$lat.utm 
All_pop_data <- subset(All_pop_data, select = -c(lon.utm, lat.utm))

All_pop_data <- transform_coord(All_pop_data, lon = "Long.x", lat = "Lat.x", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

png("./Figures_data/Plots/Admix_map5_take1_shiftedlat_long_avg.png", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon_shift.utm, y = lat_shift.utm, group = Pop, r=200000), data = All_pop_data, cols = Pop, size = 0.9) +
  geom_point(aes(x = lon_shift.utm, y = lat_shift.utm), data = All_pop_data, col="white", size=17) +
  geom_point(aes(x = lon.utm, y = lat.utm), data = All_pop_data, col="green", size=6) +
  scale_fill_manual(values=map_colours_5g)+
  geom_text(data = All_pop_data[-which(All_pop_data$Pop=="AlexOld"), ], aes(x = lon_shift.utm, y = lat_shift.utm, label = Pop), color = "black", fontface = 2, size = 25.4/72.27*20)
dev.off()

######################################





# Kluane
#check correlation with %BC and elevation

All_samples_data
Kluane <- All_samples_data[which(All_samples_data$Pop.x=="KL"|All_samples_data$Pop.x=="PC"),]

# Saximontana is V5 for Take 1

mod1 = lm(V3~Elevation, data = Kluane)
modsum = summary(mod1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.292e+00  1.440e-01   8.971 7.42e-08 ***
#   Elevation   -6.859e-04  9.243e-05  -7.421 1.00e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07808 on 17 degrees of freedom
# Multiple R-squared:  0.7641,	Adjusted R-squared:  0.7502 
# F-statistic: 55.07 on 1 and 17 DF,  p-value: 9.995e-07

png("./Figures_data/Plots/Elevation_BC_5.png", width = 3000, height = 1400)
 par(mar=c(20,20,4,4))
 plot(Kluane$Elevation, Kluane$V3, pch=20, cex = 10, mgp=c(10,5,0), cex.lab=5, cex.axis=5, xlab="Elevation(m)", ylab="BC Admixture Proportion")
 abline(mod1)
dev.off()


orderedKluane = Kluane[order(Kluane$Elevation),]

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

png("./Figures_data/Plots/Admix_Kluane_bar5.png", width = 1000, height = 707)
barplot(t(as.matrix(orderedKluane[,c(3:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(orderedKluane$ID_code), las=2, cex.names=1.7)
dev.off()


#################################################
# Climate, ice and nucleotide diversity
# Icetime.x
# tavglong
# SumWindPi

All_pop_data_rmKL.PC.GEN <- All_pop_data[-which(All_pop_data$Pop=="KL"|All_pop_data$Pop=="GEN"|All_pop_data$Pop=="PC"),]

# Models
mod3 = lm(SumWindPi ~ tavglong, data = All_pop_data_rmKL.PC.GEN)
summary(mod3)

# adding in Kluane makes it significant
#Multiple R-squared:  0.1048,	Adjusted R-squared:  0.07589 
#F-statistic: 3.628 on 1 and 31 DF,  p-value: 0.06614

# Multiple R-squared:  0.1152,	Adjusted R-squared:  0.08918 
# F-statistic: 4.427 on 1 and 34 DF,  p-value: 0.04285


mod4 = lm(Ice_time.x ~ tavglong, data = All_pop_data_rmKL.PC.GEN)
summary(mod4)

#Multiple R-squared:  0.0237,	Adjusted R-squared:  -0.00419 
#F-statistic: 0.8498 on 1 and 35 DF,  p-value: 0.3629

# Multiple R-squared:  0.01962,	Adjusted R-squared:  -0.006184 
# F-statistic: 0.7603 on 1 and 38 DF,  p-value: 0.3887


mod5 = lm(SumWindPi ~ Ice_time.x, data = All_pop_data_rmKL.PC.GEN)
summary(mod5)

#Multiple R-squared:  0.3282,	Adjusted R-squared:  0.3078 
#F-statistic: 16.12 on 1 and 33 DF,  p-value: 0.0003224

# Multiple R-squared:  0.1463,	Adjusted R-squared:  0.1226 
# F-statistic: 6.171 on 1 and 36 DF,  p-value: 0.01777

#---------------------------------
# Plots
png("./Figures_data/Plots/WindPivsTemp_sub.png", width = 1000, height = 707)
plot(All_pop_data_rmKL.PC.GEN$SumWindPi ~ All_pop_data_rmKL.PC.GEN$tavglong, col="lightblue", pch=19, cex=2, xlab="Average Temperature (degrees Celsius)", ylab="Sites Pi")
abline(mod3, col="red", lwd=3)
text(SumWindPi ~ tavglong, labels=Pop, data=All_pop_data_rmKL.PC.GEN, cex=0.9, font=2)
dev.off()

png("./Figures_data/Plots/Ice_temp_sub.png", width = 1000, height = 707)
plot(All_pop_data_rmKL.PC.GEN$Ice_time.x ~ All_pop_data_rmKL.PC.GEN$tavglong, col="lightblue", pch=19, cex=2, xlab="Average Temperature (degrees Celsius)", ylab="Icetime")
abline(mod4, col="red", lwd=3)
text(Ice_time.x ~ tavglong, labels=Pop,data=All_pop_data_rmKL.PC.GEN, cex=0.9, font=2)
dev.off()

png("./Figures_data/Plots/IcetimesWIndPi_sub.png", width = 1000, height = 707)
plot(All_pop_data_rmKL.PC.GEN$SumWindPi ~ All_pop_data_rmKL.PC.GEN$Ice_time.x, col="lightblue", pch=19, cex=2, xlab="Ice retreat time (thousands of years)", ylab="Sites Pi")
abline(mod5, col="red", lwd=3)
text(SumWindPi ~ Ice_time.x, labels=Pop, data=All_pop_data_rmKL.PC.GEN, cex=0.9, font=2)
dev.off()








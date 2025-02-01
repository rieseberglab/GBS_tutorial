# Genotyping by sequencing tutorial

authors: Cassandra Elphinstone and Yue Yu \
level: Beginner \
title: GBS tutorial \
date: Jan 2025 

# Background

This is a tutorial to walk you through how to analyze GBS data. We will cover  logging into the Digital Research Alliance servers, demultiplexing, references, mapping, SNP calling and plotting your data.

----

# Programs used in this tutorial

Below are a list of programs that we will use in this tutorial. 

### Perl 

Perl scripts will be used to demultiplex the data

### FastQC

Quality check your data


### dDocent 

dDocent is a pipeline that will walk us through building a de novo reference, mapping our reads to that reference and calling SNPs off that reference. 
It uses Rainbow to build the reference,  and FreeBayes or GATK for SNP calling.

### Admixture

Admixture is a program to look at population structure. It is a good way to tell if your data has been demultiplexed correctly. 

### R

R is a programming language that we will use to run a principal component analysis and plot our data.
Libraries needed include: ggplot2, 


----

# Useful Links

- Population genomics analyses: [https://github.com/celphin/Population_genomics_Cassiope](https://github.com/celphin/Population_genomics_Cassiope)

- Genomic analysis: https://www.zoology.ubc.ca/~schluter/R/Genomics.html#SNP-calling_pipeline



----

# Zoom Recordings of past workshops
2025 January 30th and 31st:
https://drive.google.com/drive/folders/1irEqBinzt6u-6G7kAR2fvnVIeNml3aQD?usp=sharing


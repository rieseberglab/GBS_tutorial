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
[Jan 30th recording - Part 1:](
https://drive.google.com/drive/folders/1irEqBinzt6u-6G7kAR2fvnVIeNml3aQD?usp=share_link)

[Jan 31st recording - Part 2.1:](
https://ubc.zoom.us/rec/share/9jz05It4nh05hhfy4j8-O5_KfQfX1Hfu3meamg74wjDMXItsrb-LZ4unTXL_4FNq.wmBRbsBp11CgZH7A)
Password: G?n$@=9Q

[Jan 31st recording - Part 2.2:](
https://ubc.zoom.us/rec/share/Lfyy9fVSBUL4HfeBRYS4QXB3-vRMbF6P-e6LeLhZPNF87bo4ld-sM_7xSA1bVHx3.UuYcJ6UZB5xqItTo)
Password:*7.yVqU4

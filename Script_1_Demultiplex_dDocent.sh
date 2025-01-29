#####################################
#
# Part 1: GBS data demultiplex & dDocent
# January 2025
#
# Code Contributor: Cassandra E.
# Biodiversity Research Center, UBC
#
#####################################


# This code includes the following 3 sections:
# 1.1: Demultiplexing
# 1.2: Use of dDocent (Mapping -> SNP calling -> SNP filtering)
# 1.3: Example code on how to run dDocent


#####################################
#  1.1 Demultiplexing 
#####################################

# -- Move the following 3 files into your working directory
mkdir working_dir

# -- File 1: fastq file (R1 and R2)
cp sample.R1.fastq.gz working_dir
cp sample.R2.fastq.gz working_dir

# -- File 2: barcode file
cp barcodes.txt working_dir

# -- File 3: perl script
cp GBS2enzymedemultiplex_notrim.pl working_dir





# -- TMUX: allows to run in the background
tmux new-session -s lane1
tmux attach-session -t lane1

# -- Start interactive sessions (SALLOC) to run this step
cd working_dir
salloc -c2 --time 20:00:00 --mem 120000m --account def-rieseber

# -- Demultiplexing
perl GBS2enzymedemultiplex_notrim.pl barcodes.txt sample.R1.fastq.gz sample.R1.fastq.gz 






#####################################
#  1.2 dDocent
#####################################

# -- create tmux session 
tmux new-session -s Cassiope
tmux attach-session -t Cassiope

# -- get interactive allocation
salloc -c48 --time 23:00:00 --mem 187G --account def-rieseber

# -- set up Ddocent environment in session
~/miniconda2/bin/conda create --name dDocent  
source ~/miniconda2/bin/activate dDocent

# -- update conda
conda update -n base -c defaults conda
# Your installed version is: 2.17

#source ~/miniconda2/bin/deactivate dDocent
conda install -c bioconda ddocent
# Update from dDocent 2.6.0 to dDocent 2.9.4




##########################################
#  1.3 Example code on how to run dDocent
##########################################

#----------------
# change directory to your fastq files named according to Ddocent:
# https://www.ddocent.com/UserGuide/#quality-filtering

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/reference_fastq

# activate conda environment
source ~/miniconda2/bin/activate dDocent

#use dDocent to build the assembly (leave k2 parameter blank)
dDocent

55 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes

dDocent detects 48 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
48

dDocent detects 187 gigabytes of maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes
For example, to limit dDocent to ten gigabytes, enter 10
This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0
Then press [ENTER]
180

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
yes

What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap.
Then press [ENTER]
PE

Reads will be assembled with Rainbow

CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa.
Would you like to enter a new c parameter now? Type yes or no and press [ENTER]
yes
Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)
0.95

Do you want to map reads?  Type yes or no and press [ENTER]
no

Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
cassandra.elphinstone@shaw.ca

dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background.


                       Number of Unique Sequences with More than X depth (Counted within individuals)

    3e+06 ++----------+-----------+----------+-----------+-----------+-----------+----------+-----------+----------++
          +           +           +          +           +           +           +          +           +           +
          *                                                                                                         |
          |*                                                                                                        |
  2.5e+06 ++*                                                                                                      ++
          | *                                                                                                       |
          |  *                                                                                                      |
          |   *                                                                                                     |
    2e+06 ++   *                                                                                                   ++
          |    *                                                                                                    |
          |     *****                                                                                               |
  1.5e+06 ++         *                                                                                             ++
          |           *****                                                                                         |
          |                *                                                                                        |
          |                 ******                                                                                  |
    1e+06 ++                      *****                                                                            ++
          |                            ************                                                                 |
          |                                        ******************                                               |
          |                                                          ******************                             |
   500000 ++                                                                           *****************************+
          |                                                                                                         *
          |                                                                                                         |
          +           +           +          +           +           +           +          +           +           +
        0 ++----------+-----------+----------+-----------+-----------+-----------+----------+-----------+----------++
          2           4           6          8           10          12          14         16          18          20
                                                           depth

Please choose data cutoff.  In essence, you are picking a minimum (within individual) coverage level for a read (allele) to be used in the reference assembly
3



                                 Number of Unique Sequences present in more than X Individuals
 Number of Unique Sequences
   300000 ++----------------+----------------+-----------------+-----------------+----------------+----------------++
          +                 +                +                 +                 +                +                 +
          |      *                                                                                                  |
          |       *                                                                                                 |
   250000 ++      *                                                                                                ++
          |        *                                                                                                |
          |        *                                                                                                |
          |         *                                                                                               |
   200000 ++        *                                                                                              ++
          |          *                                                                                              |
          |           *                                                                                             |
   150000 ++           *                                                                                           ++
          |            *                                                                                            |
          |             ***                                                                                         |
          |                *                                                                                        |
   100000 ++                ***                                                                                    ++
          |                                                                                                         |
          |                    ****                                                                                 |
          |                        ***                                                                              |
    50000 ++                          ****                                                                         ++
          |                               ***                                                                       |
          |                                  ****                                                                   |
          +                 +                +   ******************              +                +                 +
        0 ++----------------+----------------+-----------------+---***************************************---------++
          0                 5                10                15                20               25                30
                                                     Number of Individuals

Please choose data cutoff.  Pick point right before the assymptote. A good starting cutoff might be 10% of the total number of individuals
5


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type 'bg' without the quotes and press enter
Type 'disown -h' again without the quotes and press enter

Now sit back, relax, and wait for your analysis to finish

dDocent assembled 102787 sequences (after cutoffs) into 23890 contigs

dDocent has finished with an analysis in /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq/reference_fastq

dDocent started Mon Mar 6 16:58:48 PST 2023

dDocent finished Mon Mar 6 17:04:39 PST 2023

dDocent 2.9.4
The 'd' is silent, hillbilly.



# END

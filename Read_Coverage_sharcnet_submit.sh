#!/bin/bash

# SHARCNet submission script to get read coverage per individual per locus.

# Requires:
# 1 - vcf file in working directory. Ideally filtered for missing individuals
# missing loci, and minor allele frequency.
# Output format has to be identical to filename after --out (see form below)
# vcftools --vcf  variants_rbtassem_0.5_maf0.05.recode.vcf\
# --extract-FORMAT-info AD --out variants_gbct.recode
# 2 - partial version of the filter paraolgs script:
# partial_filter_paralogs.R

### ---------- Job configuration --------------------------------------------

# Run dependent and permanent parameters
# will be run on complete nodes NOT partial

#SBATCH --nodes=1                       # number of nodes to use
#SBATCH --time=00-00:05:00              # time (DD-HH:MM:SS)
#SBATCH --account=def-emandevi          # account name
#SBATCH --job-name="Read coverage - miss0.5 maf0.05 S" # name of job

#SBATCH --ntasks-per-node=32            # taks per node (one core per node)
#SBATCH --mem=128000M                   # memory per node
#SBATCH --output=RC[0.5_0.05_S].log             # log file
#SBATCH --error=RC[0.5_0.05_S].err              # error file
#SBATCH --mail-user=bratz@uoguelph.ca   # who to email
#SBATCH --mail-type=END                 # when to email

# Modules to Load:
# Only need R and prerequests
module load nixpkgs/16.09 gcc/7.3.0
module load r/3.5.1

### ---------- Functions ----------------------------------------------------

# USE: function Dead() kills the process if files are not present in your
# path or if something goes wrong
dead () {

  # Prints warning message that process is about to be terminated
  echo "********** WARNING: EARLY TERMINATION **********"
  echo "Something went wrong...oops"

  # Quits running program
  exit 1
}

### --------- Main ----------------------------------------------------------

# Get working directory
WORKDIR=$PWD

# Loop to check if request files and Rscript are in your working directory
if test -f  variants_gbct.recode.AD.FORMAT || test -f partial_filter_paralogs.R; then
  # File was found.
  # Assigning to a variable
  IN_FILE="${WORKDIR}/variants_gbct.recode.AD.FORMAT"

  # Script was found
  # Assigning to a variable
  OUTSIDE_SCRIPT="${WORKDIR}/partial_filter_paralogs.R"

else
  # File wasn't found.
  # Dead function called to terminate the program early
  echo -en "missing something in your path.\n"
  dead

fi

# Removing commas between alleles using sed
echo "Removing commas between alleles"
sed 's/\/project\/ysctrout\/emandevi\/gbcutthroat\/bwa_assem_rbt\/aln_//g' $IN_FILE | sed 's/\.sorted\.bam//g' | sed 's/\t/,/g' > variants_simplifiedIDs.txt

# Echo starting to calculate read count using R script
echo """
  Calling R script.
  Calculating reads per individual per locus.
"""

# Call R script
Rscript $OUTSIDE_SCRIPT

# Echo finished calculating read count
echo """
  Finished calculating reads per individual per locus.
  Written to file 'readsperindperlocus.txt'.
  ALso echoed below to final count also appears in log file.
"""

# loop that finds output file and echoes the number of reads per locus
# so the result also appears in the log file.
if test -f readsperindperlocus.txt; then
  # Set file to variable
  ANSWER="${WORKDIR}/readsperindperlocus.txt"

  # echo the answer for the log file
  echo $ANSWER

else
  # Wasn't able to find the file. Something went wrong
  # Dead fucntion calles to terminate the program
  echo -en "no file could be found.\nTry again."
  dead

fi

# Echo program done
echo "Program done. Goodbye"

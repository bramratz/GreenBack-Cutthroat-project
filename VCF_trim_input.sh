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

#SBATCH --nodes=1                                           # number of nodes to use
#SBATCH --time=00-03:00:00                         # time (DD-HH:MM:SS)
#SBATCH --account=def-emandevi               # account name
#SBATCH --job-name="GBCT vcf filtering" # name of job

#SBATCH --ntasks-per-node=32            # taks per node (one core per node)
#SBATCH --mem=128000M                   # memory per node
#SBATCH --output=sim-%j.log             # log file
#SBATCH --error=sim-%j.err              # error file
#SBATCH --mail-user=bratz@uoguelph.ca   # who to email
#SBATCH --mail-type=END                 # when to email

# Load required modules - vcftools
# ***need to load R with version attached - get error otherwise***
module load nixpkgs/16.09 intel/2016.4
module load nixpkgs/16.09 intel/2018.3
module load vcftools/0.1.14

### ---------- Command line job arguments -----------------------------------
# Outlines command line arguments that are required in order to set up
# vcftools pipe. If arguments are not specified a usage function will be
# displayed and the program will exit back to the commandline prompt
# Usage function outlining what arguments to enter. If nothing entered will
# be displayed
usage () {
  echo -e "\nUsage: $0 [-v <vcf file>][-l <missing data per locus>] [-i <missing data per individual>] [-f <minor allele frequency>] [-p <paralog filtering>]" 1>&2;
  echo -e """
  Program arguments:
    -v Requires a vcf file
    -l Requires a float between 0 and 1 specifying the minimum amount of data
    a locus must have for the SNP to be identified.
    -i Requires a float between 0 and 1 specifying the maximum amount of loci
    an individual can be missing to be retained.
    -f Requires a float between 0 and 1 specifying sites to be retained only if
    they have a minor allele frequency greater than or equal to the value
    given. If no maf filter desired set to 0.
    -p Requires a one of three string options: none, permissive, strict,
  """
  exit 1;
}

# Handles the command line inputs. Assigns each input to a variable, if no
# variable provided will echo which variable is missing and call the usage
# usage function.
while getopts ":v:l:i:f:p:" o; do
  case "${o}" in
    v )
      v=${OPTARG}
      if [[ -z "${v}" ]]
        then
          echo -en "Missing vcf file"
          usage
      fi
      ;;
    l )
      l=${OPTARG}
      if [[ -z "${l}" ]]
        then
          echo -en "Missing value for l"
          usage
      fi
      ;;
    i )
      i=${OPTARG}
      if [[ -z "${i}" ]]
        then
          echo -en "Missing value for i"
          usage
      fi
      ;;
    f )
      f=${OPTARG}
      if [[ -z "${f}" ]]
        then
          echo -en "Missing value for f"
          usage
      fi
      ;;
    p )
      p=${OPTARG}
      if [[ -z "${p}" ]]
        then
          echo -en "Missing string for p"
          usage
      fi
      ;;
#    * )
#      usage
#      ;;
  esac
done
shift $((OPTIND-1))

# Further check to ensure the arguments were inputed from the command line
# ***Possibily redundant***
if [ -v "${v}" ] || [ -z "${l}" ] || [ -z "${i}" ] || [ -z "${f}" ]; then
  usage
fi

### ---------- Useful job infoformation -------------------------------------

echo "Current working directory: $PWD"
echo "starting run at: $DATE"
echo "------------------------------------------------"
echo "job is running on node: $HOSTNAME"
echo "------------------------------------------------"
echo "Job identifier is $SLURM_JOB_ID"
echo "Job name is $SLURM_JOB_NAME"

# Echo variables to make sure they are correctly assigned
echo "v = ${v}"
echo "l = ${l}"
echo "i = ${i}"
echo "f = ${f}"
echo "p = ${p}"


### ---------- Filtering options --------------------------------------------

# USE: function Dead kills the process if files have not been created that are
# required for further filtering
dead () {

  # Prints warning message that process is about to be terminated
  echo "********** WARNING: EARLY TERMINATION **********"
  echo "File not found, terminating process."

  # Quits running program
  exit 1
}

# USE: fucntion filter_raw uses command line inputs as values for arguments
# in the vcftools function. The function pipes together multiple vcftools
# calls in order to filter everything specified from the command line with the
# exception of paralogs (handeled separately).
# Filters: missing data per locus:0.5 ; missing data for individual:>95% ; minor allele freq:0.05
filter_raw () {

  # First filtering step - only keep variants that have been successfully
  # genotyped in > 50% of individuals AND have a minor allele frequency greater
  # than 0.05
  vcftools --vcf $v --max-missing $l --maf $f --recode --recode-INFO-all --out raw_miss05_maf005

  # Second step - Remove individuals that did not sequence well, if individual
  # has > 95% data missing it is removed
  # Create file of missing data for each individual
  vcftools --vcf raw_miss05_maf005.recode.vcf --missing-indv

  # Filter out individual if missing > 95%
  awk '$5 > 0.95' out.imiss | cut -f1 > rm_indv_gt95missing_maf005

  # Filter out individuals in list from raw vcf file
  vcftools --vcf raw_miss05_maf005.recode.vcf --remove rm_indv_gt95missing_maf005 --recode --recode-INFO-all --out filtered_pre_paralog
}

### ---------- Main ---------------------------------------------------------

# Assign path for raw vcf file to a varaible
#RawData="/home/bramratz/projects/def-emandevi/brmaratz/RaW/gbcutthroat_rbtassem_rawvariants.vcf"
paralog_script='/home/bramratz/projects/def-emandevi/bramratz/filter_paralogs.R'
perl_script1='/home/bramratz/projects/def-emandevi/bramratz/vcf2mpgl.pl'
perl_script2='./gl2genest.pl'
workdir=$(pwd)

# Call filter_raw function to begin filtering vcf file. Output saved as
# filtered_pre_paralog variable
#filtered_pre_paralog=$(filter_raw $RawData $l $i $f)

# test function
filter_raw $v $l $i $f

# House cleaning
# use file test operators -f,-s
if test -f filtered_pre_paralog.*; then

  # Saves file name as a variable
  Good_file=filtered_pre_paralog.*

  # Prints file name
  echo "file" $Good_file "exists and is of non-zero size."
  echo "Removing intermediate files."

  # Compresses intermediate files with gzip and moves to separate directory
  # if needed later
  gzip raw_miss05_maf005.recode.vcf out.imiss rm_indv_gt95missing_maf005

  mv -t ./inter_files raw_miss05_maf005.recode.vcf.gz\
    out.imiss.gz\
    rm_indv_gt95missing_maf005.gz

else
  # File not found message
  echo "Required file does not exist or is of file size zero."

  # Terminates running process
  dead

fi

# Convert vcf file into file format that filter_paralogs program can use
vcftools --vcf filtered_pre_paralog.recode.vcf --extract-FORMAT-info AD --out variants_gbct.recode

# Remove the commas between alleles with sed
sed 's/\/project\/ysctrout\/emandevi\/gbcutthroat\/bwa_assem_rbt\/aln_//g' variants_gbct.recode.AD.FORMAT | sed 's/\.sorted\.bam//g' | sed 's/\t/,/g' > variants_miss0.5_maf0.05_ind0.95alleledepth_simplifiedIDs.txt

# load perl module
module load perl/5.22.4

# Generate mgpl file for paralog_filtering script using vcf2mpgl.pl
perl $perl_script1 filtered_pre_paralog.recode.vcf

# Add headers
cat header_mpglfiltered_pre_paralog.recode.txt filtered_pre_paralog.recode.mpgl > gl2genes_pre_paralog.mpgl

# ***altered the gl2genest.pl file in order to set the arguments up correctly
# the code was taking the script is 1 instead of zero for some reason

# Take that file and put it into another perl program - I have no idea if [1/1] is right but that's what I used
## fix this via hard coding
#perl $perl_script2 gl2genes_pre_paralog.mpgl mean
perl $perl_script2 ./gl2genes_pre_paralog.mpgl mean
# output is "pntest_"."$summarytype"."_$out"
# *** will have to change R script to handle new output file

# Load modules
module load nixpkgs/16.09 gcc/7.3.0
module load r/3.5.1

# Read allele depth file into filter paralogs R script
Rscript $paralog_script

# Reload vcftools
module load nixpkgs/16.09 intel/2016.4
module load nixpkgs/16.09 intel/2018.3
module load imkl/11.3.4.258
module load vcftools/0.1.14

# For loop for determining which of the two parolog filtering options
# (permissive or strict) should be used based on command line arguments
if [ $p = "strict" ]; then
  rm keep_pos_permissive.txt
  vcftools --vcf filtered_pre_paralog.recode.vcf --remove keep_pos_strict.txt --recode --recode-INFO-all --out final_filtered_strict
elif [ $p == "permissive" ]; then
  rm keep_pos_strict.txt
  vcftools --vcf filtered_pre_paralog.recode.vcf --remove keep_pos_permissive.txt --recode --recode-INFO-all --out final_filtered_permissive
elif [ $p == "none" ]; then
  pass
else
  dead
fi

# recall vcf2mpgl on final filter to prepare for entropy
# loop that removes previous files so names don't get confusing
if test -f final_filtered_*; then
  # Compresses intermediate files with gzip and moves to separate directory
  # if needed later
  gzip filtered_pre_paralog.recode.vcf\
    pntest_mean_gl2genes_pre_paralogs.txt\
    header_mpglfiltered_pre_paralog.recode.txt\
    filtered_pre_paralog.mpgl\
    gl2genes_pre_paralog.mpgl\
    indiv_ids.txt

  mv -t ./inter_files filtered_pre_paralog.recode.vcf.gzip\
    pntest_mean_gl2genes_pre_paralogs.txt.gzip\
    header_mpglfiltered_pre_paralog.recode.txt.gzip\
    filtered_pre_paralog.mpgl.gzip\
    gl2genes_pre_paralog.mpgl.gzip\
    indiv_ids.txt.gzip

   # use vcf2mpgl.pl again
   perl $perl_script1 final_filtered_strict.recode.vcf

   # Add headers
   cat header_mpglfinal_filtered_strict.recode.txt final_filtered_strict.recode.mpgl > entropy_in_variants_all.mpgl

else
  # files failed to be created kill process
  dead

fi

# end of filtering

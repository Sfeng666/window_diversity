#!/bin/bash

### 1. set paths ###
## set parameter parsings for user-specified paths ##
usage="## this script calculate window nucleotide diversity.\n
usage: bash $0 [-h/--help] [-w/--work_directory <path>] [-s/--dir_snp <path>]\n
[-m/--dir_mpileup <path>] [-c/--chr_assignment <path>] [-a/--sample_size <path>]\n
[-g/--window_heter <int>] [-v/--mincov <int>] [-x/--maxcov <path>] [-t/--min_count <int>]\n
[-l/--calc_script <path>] [-p/--script <path>] [-d/--conda_path <path>] [-e/--conda_env <name>]"

while [ $# -gt 0 ]; do
  case "$1" in
    --work_directory*|-w*)
      if [[ "$1" != *=* ]]; then shift; fi # Value is next arg if no `=`
      wd="${1#*=}"
      ;;
    --dir_snp*|-s*)
      if [[ "$1" != *=* ]]; then shift; fi
      dir_snp="${1#*=}"
      ;;
    --dir_mpileup*|-m*)
      if [[ "$1" != *=* ]]; then shift; fi
      dir_mpileup="${1#*=}"
      ;;
    --chr_assignment*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      chr_assignment="${1#*=}"
      ;;
    --sample_size*|-a*)
      if [[ "$1" != *=* ]]; then shift; fi
      sample_size="${1#*=}"
      ;;
    --window_heter*|-g*)
      if [[ "$1" != *=* ]]; then shift; fi
      window_heter="${1#*=}"
      ;;
    --min_cov*|-v*)
      if [[ "$1" != *=* ]]; then shift; fi
      min_cov="${1#*=}"
      ;;
    --max_cov*|-x*)
      if [[ "$1" != *=* ]]; then shift; fi
      max_cov="${1#*=}"
      ;;
    --min_count*|-t*)
      if [[ "$1" != *=* ]]; then shift; fi
      min_count="${1#*=}"
      ;;
    --calc_script*|-l*)
      if [[ "$1" != *=* ]]; then shift; fi
      calc_script="${1#*=}"
      ;;
    --script*|-p*)
      if [[ "$1" != *=* ]]; then shift; fi
      script="${1#*=}"
      ;;
    --conda_path*|-d*)
      if [[ "$1" != *=* ]]; then shift; fi
      conda_path="${1#*=}"
      ;;
    --conda_env*|-e*)
      if [[ "$1" != *=* ]]; then shift; fi
      conda_env="${1#*=}"
      ;;      
    --help|-h)
      echo -e $usage # Flag argument
      exit 0
      ;;
    *)
      >&2 printf "Error: Invalid argument\n"
      exit 1
      ;;
  esac
  shift
done

if [[ -z $wd || -z $dir_snp || -z $dir_mpileup || -z $chr_assignment  \
|| -z $sample_size || -z $window_heter \
|| -z $max_cov || -z $min_cov || -z $min_count || -z $calc_script || -z $script \
|| -z $conda_path || -z $conda_env ]]; then
    echo -e $usage
    echo "-w $wd
-s $dir_snp
-m $dir_mpileup
-c $chr_assignment
-a $sample_size
-g $window_heter
-v $min_cov
-x $max_cov
-t $min_count
-l $calc_script
-p $script
-d $conda_path
-e $conda_env"
    exit 1
fi

## derived paths ##
log=$wd/calc_Fst.log
nt_chr=$wd/diversity_{2}.txt
window_stats=$wd/window_stats.txt
temp=$wd/temps
snp_clean_contig=$dir_snp/{1}_snp_clean.vcf
indel_positions_contig=$dir_snp/{1}_inDel-positions_20.txt
mpileup_contig=$dir_mpileup/{1}_mpileup.gz

mkdir -p $wd
mkdir -p $temp

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" 'Chr' 'Contig' 'Window' 'Start' 'End' 'Length' 'Number of filtered sites' 'Number of SNPs' > $window_stats

names=$(awk 'BEGIN{FS = "\t"; OFS = ","} NR == 1{print $0}' $sample_size)

### 2. configure the conda enviroment ###
source $conda_path/bin/activate $conda_env

### 3. run steps of the pipeline ###
### job starts ###
echo start at time `date +%F'  '%H:%M`

echo -e "### Whole job starts ###\n" > $log

## define a function to calculate window nucleotide diversity for each contig
calc_Fst_contig() {
local contig=$1
local chr=$2
local snp_clean_contig=$3
local indel_positions_contig=$4
local mpileup_contig=$5
local nt_chr=$6

python $calc_script \
--in_snp $snp_clean_contig \
--in_mpileup $mpileup_contig \
--in_indel_pos $indel_positions_contig \
--in_sample_size $sample_size \
--out_window_stats $window_stats \
--out_window_diversity $nt_chr \
--script $script \
--window_heter $window_heter \
--min_count $min_count \
--min_cov $min_cov \
--max_cov $max_cov \
--contig $contig \
--chr $chr \
--names $names \
--temp $temp

echo "nucleotide diversity calculated from $snp_clean_contig at chromosome arm $chr" >> $log
}

## For each contig,

## 3.1. calculate the window nucleotide diversity for each contig.

export -f calc_Fst_contig
export sample_size
export window_stats
export calc_script
export script
export window_heter
export min_count
export min_cov
export max_cov
export names
export temp
export log

echo -e "### Estimation of window-based nucleotide diversity starts ###\n" >> $log

awk 'BEGIN{OFS="\t"} $2 != "na" && NR>1{print $1, $2}' $chr_assignment | sort -k1 \
| parallel --tmpdir $temp -j 0 -k --colsep '\t' \
calc_Fst_contig {1} {2} $snp_clean_contig $indel_positions_contig $mpileup_contig $nt_chr \

echo -e "### Estimation of nucleotide diversity ends ###\n" >> $log
echo -e "### Whole job done ###\n" >> $log
rm -r $temp

### script finishes ###
echo finish at time `date +%F'  '%H:%M`
from optparse import OptionParser
import os
from re import T

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: calculate window nucleotide diversity'''
version = '%prog 09.05.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--in_snp",
                    action="store",
                    dest = "in_snp",
                    help = "input path of filtered SNPs from a specified chr in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--in_mpileup",
                    action="store",
                    dest = "in_mpileup",
                    help = "input path of mpileup file of a specified chr in mpileup format (.mpileup)",
                    metavar = "PATH")
parser.add_option("--in_indel_pos",
                    action="store",
                    dest = "in_indel_pos",
                    help = "input path of the file containing indel positions (.txt)",
                    metavar = "PATH")
parser.add_option("--in_sample_size",
                    action="store",
                    dest = "in_sample_size",
                    help = "input path of allelic (auto/X) sample sizes of all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--out_window_stats",
                    action="store",
                    dest = "out_window_stats",
                    help = "output path of window stats of all chrs for all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--out_window_diversity",
                    action="store",
                    dest = "out_window_diversity",
                    help = "output path of window nucleotide diversity of each chromosome arm (.txt)",
                    metavar = "PATH")
parser.add_option("--script",
                    action="store",
                    dest = "script",
                    help = "path of the python script counting the number of sites within each window (.py)",
                    metavar = "PATH")
parser.add_option("--window_heter",
                    action="store",
                    dest = "window_heter",
                    help = "the threshold of accumulative average heterozygosity among all samples analyzed, that is used to define windows",
                    metavar = "INT")
parser.add_option("--min_count",
                    action="store",
                    dest = "min_count",
                    help = "the threshold of minimum minor allele count (to exclude low-count alleles that might be sequencing errors)",
                    metavar = "FLOAT")                    
parser.add_option("--min_cov",
                    action="store",
                    dest = "min_cov",
                    help = "the threshold of minimum coverage",
                    metavar = "FLOAT")
parser.add_option("--max_cov",
                    action="store",
                    dest = "max_cov",
                    help = "input file containing the maximum coverage threshold of each contig",
                    metavar = "PATH")
parser.add_option("--contig",
                    action="store",
                    dest = "contig",
                    help = "a string indicating the contig being processed",
                    metavar = "STR")
parser.add_option("--chr",
                    action="store",
                    dest = "chr",
                    help = "a string indicating which chromosome arm in D.mel the contig uniquely maps to",
                    metavar = "STR")
parser.add_option("--names",
                    action="store",
                    dest = "names",
                    help = "a comma-delimited string of sample names",
                    metavar = "STR")
parser.add_option("--temp",
                    action = "store",
                    dest = "temp",
                    help = "output directory of temporary intermediate files",
                    metavar = "PATH")
(options,args) = parser.parse_args()

### introduced variables ###
in_snp = options.in_snp
in_mpileup = options.in_mpileup
in_indel_pos = options.in_indel_pos
in_sample_size = options.in_sample_size
out_window_stats = options.out_window_stats
out_window_diversity = options.out_window_diversity
script = options.script
window_heter = int(options.window_heter)
min_cov = options.min_cov
max_cov = options.max_cov
min_count = options.min_count
contig = options.contig
chr = options.chr
names = options.names
temp = options.temp

### def functions ###

## function to calculate single SNP heterozygosity
def ht_SNP(ct):
    # calcualte from the two alleles that have the top-two allele count
    ct_top2 = sorted(ct, reverse = True)[:2]
    dp = sum(ct_top2)
    ht = (1 - (ct_top2[0]/dp)**2 - (ct_top2[1]/dp)**2)*dp/(dp - 1)
    return ht

### 1. read the input of sample size
sample_size = {'auto': [], 'X': []}
with open(in_sample_size, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i == 0:
            samples = line
        elif i == 1:
            sample_size['auto'] = {x: int(y) for x,y in zip(samples, line)}
        elif i == 2:
            sample_size['X'] = {x: int(y) for x,y in zip(samples, line)}
        i += 1
if chr == 'X':
    sample_size_contig = sample_size['X']
else:
    sample_size_contig = sample_size['auto']

### 2. partition windows based on summed heterozygosity and calculate pairwise Fst within windows
windows = {}
win, sum_heter, num_snp = 0, 0, 0
windows[win] = {'start': 1, 
'accum_ht': {sample: 0 for sample in samples}}

start = 0
end = 0
with open(in_snp, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            num_snp += 1
            pos = int(line[1])
            ref = line[3].upper()
            alt = line[4].upper()
            alleles = [ref] + alt.split(',')
            afs = []
            dps = []

            # specify the ending boundry of each window at the midpoint of the end SNP and the next starting SNP
            if end == 1 and start == 1:
                sum_pos = windows[win]['end'] + pos
                if sum_pos & 1 == 0:
                    windows[win]['end'] = int(sum_pos/2)
                else:
                    windows[win]['end'] = int((sum_pos - 1)/2)
                win += 1
                windows.setdefault(win, {'start': windows[win - 1]['end'] + 1,
                'accum_ht': {sample: 0 for sample in samples}})
                end = 0
                start = 0

            # specify the starting boundry of each window at the midpoint of the starting SNP and the previous ending SNP + 1
            if start == 1:
                windows.setdefault(win, {'start': windows[win - 1]['end'] + 1,
                'accum_ht': {sample: 0 for sample in samples}})
                start = 0
            
            # extract allele frequency of all alleles for all samples and calculate heterozygosity
            for i in range(len(line[9:])):
                sample = samples[i]
                info = line[9:][i]
                ct = [int(info.split(':')[1])] + list(int(x) for x in info.split(':')[2].split(','))
                windows[win]['accum_ht'][sample] += ht_SNP(ct)
                sum_heter += ht_SNP(ct)

            # partition windows based on accumulative heterozygosity vs threshold
            if sum_heter/len(samples) >= window_heter:
                windows[win]['end'] = pos
                windows[win]['num_snp'] = num_snp
                end = 1
                start = 1

                sum_heter = 0
                num_snp = 0

    if sum_heter/len(samples) < window_heter:
        if sum_heter/len(samples) >= window_heter*0.95:
            windows[win]['end'] = pos
            windows[win]['num_snp'] = num_snp
        else:
            del windows[win]

### 3. output the window positions as the input of a script counting the number of filtered sites within each window
if len(windows) > 0:
    window_pos = contig + '.pos'
    window_ct = contig + '.ct'
    with open(window_pos, 'w') as f:
        for win in windows:
            f.write('\t'.join([str(windows[win]['start']), str(windows[win]['end'])]) + '\n')

    ### 4. run the script to count the number of filtered sites within each window
    ct_script = contig + '.sh'
    indelfree_mpileup = contig + '.mpileup'
    with open(ct_script, 'w') as f:
        f.write('''awk 'BEGIN{FS="\\t"} \
NR==FNR {D[$1$2]++; next} \
!($1$2 in D) && NR>FNR {print}' \
<(cat ''' + in_indel_pos + ''') \
<(zcat ''' + in_mpileup + ''') \
> ''' + indelfree_mpileup + '''
python2.7 ''' + script + ''' \
--mpileup ''' + indelfree_mpileup + ''' \
--min-cov ''' + min_cov + ''' \
--max-cov ''' + max_cov + ''' \
--min-count ''' + min_count + ''' \
--min-freq 0 \
--miss-frac 0.001 \
--names ''' + names + ''' \
--windows ''' + window_pos + ''' \
--windows_ct ''' + window_ct)

    os.system('bash ' + ct_script)
    os.system('rm ' + ct_script)
    os.system('rm ' + indelfree_mpileup)
    os.system('rm ' + window_pos)

    ### 5. output the stats of windows
    with open(out_window_stats, 'a') as f, open(window_ct, 'r') as f_ct_win:
        i = 0
        for line in f_ct_win:
            if i == 0:
                num_sites = line.strip()
            else:
                line = line.strip().split('\t')
                windows[int(line[0])]['num_sites'] = int(line[3])
                # header - 'Chr' 'Contig' 'Window' 'Start' 'End' 'Length' 'Number of filtered sites' 'Number of SNPs'
                f.write('\t'.join([chr, contig, line[0], line[1], line[2], str(int(line[2]) - int(line[1]) + 1), line[3], str(windows[int(line[0])]['num_snp'])]) + '\n')
            i += 1
    os.system('rm ' + window_ct)

    ### 6. Calculate nucleotide diversity  
    header = samples
    print_header = True
    if os.path.isfile(out_window_diversity):
        print_header = False
    with open(out_window_diversity, 'a') as f:
        if print_header:
            f.write('\t'.join(header) + '\n')
        for win in windows:
            nt_win = list(str(windows[win]['accum_ht'][sample]/windows[win]['num_sites']*sample_size_contig[sample]/(sample_size_contig[sample] - 1)) for sample in samples)
            f.write('\t'.join(nt_win) + '\n')
from optparse import OptionParser
import gzip
import re 

### help and usage ### 
usage = "usage: %prog [options] args"
description = '''Function: extract gene symbols of genes from windows of the lowest diversity'''
version = '%prog 09.05.2021'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--wd",
                    action="store",
                    dest = "wd",
                    help = "working directory",
                    metavar = "PATH")
parser.add_option("--in_win_stats",
                    action="store",
                    dest = "in_win_stats",
                    help = "input file of window stats of all chrs for all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--in_win_lowdiv",
                    action="store",
                    dest = "in_win_lowdiv",
                    help = "input file of windows that have the lowest nucleotide diversity within each window (.txt)",
                    metavar = "PATH")
parser.add_option("--in_anno",
                    action="store",
                    dest = "in_anno",
                    help = "input path of genome annotation in compressed gff3 format (.gff.gz)",
                    metavar = "PATH")
parser.add_option("--in_ortho",
                    action="store",
                    dest = "in_ortho",
                    help = "input path of suzukii gene table orthlogous to melanogaster (.txt)",
                    metavar = "PATH")
parser.add_option("--out_gene",
                    action="store",
                    dest = "out_gene",
                    help = "output file of extracted gene symbols (.txt)",
                    metavar = "PATH")
(options,args) = parser.parse_args()

### introduced variables ###
wd = options.wd
in_win_stats = options.in_win_stats
in_win_lowdiv = options.in_win_lowdiv
in_anno = options.in_anno
in_ortho = options.in_ortho
out_gene = options.out_gene

### def functions ###

### 1. read the annotation file from RefSeq gff
annotation = {}
with gzip.open(in_anno, 'rt') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            chr = line[0]
            feature = line[2]
            start = int(line[3])
            end = int(line[4])
            strand = line[6]
            attr = line[8]
            if feature == 'region':
                annotation.setdefault(chr, {})
            elif feature == 'gene':
                gene = re.search('Dbxref=GeneID:(.*?);', attr).group(1)
                try:
                    gene_biotype = re.search('gene_biotype=(.*?);', attr).group(1)
                except:
                    gene_biotype = re.search('gene_biotype=(.*)', attr).group(1)
                if gene_biotype == 'protein_coding':
                    annotation[chr].setdefault(gene, [start, end, strand])

### 2. read the orthology table
ortho = {}
with open(in_ortho, 'r') as f:
    for line in f:
            line = line.strip().split('\t')
            gene = line[0]
            symb = line[2]
            fb = line[4].split(':')[1]
            des = line[5]
            ortho[gene] = [symb, fb, des]

### 3. build a dictionary for accessions of analyzed contigs and windows
win_stats = {}
with open(in_win_stats, 'r') as f:
    i = 0
    for line in f:
        if i != 0:
            line = line.strip().split('\t')
            chr = line[0]
            contig = line[1]
            win = line[3]
            start = int(line[4])
            end = int(line[5])
            if not chr in win_stats:
                win_stats[chr] = {contig: {win: [start, end]}}
            else:
                if not contig in win_stats[chr]:
                    win_stats[chr][contig] = {win: [start, end]}
                else:
                    win_stats[chr][contig][win] = [start, end]
        i += 1

### 4. extract windows of the lowest nucleotide diversity within each contig
win_div = {}
with open(in_win_lowdiv, 'r') as f, open(out_gene, 'w') as out:
    header = ['Chromosomal arm', 'Contig', 'Window', 'Gene ID', 'Gene Symbol', 'Flybase ID', 'Description']
    out.write('\t'.join(header) + '\n')
    for line in f:
            line = line.strip().split('\t')
            contig = line[0]
            chr = list(chr for chr in win_stats if contig in win_stats[chr])[0]
            win = line[2]
            for gene in annotation[contig]:
                if annotation[contig][gene][0] >= win_stats[chr][contig][win][0] and annotation[contig][gene][1] <= win_stats[chr][contig][win][1]:
                    if gene in ortho:
                        wline = [chr, contig, win, gene, ortho[gene][0], ortho[gene][1], ortho[gene][2]]
                        out.write('\t'.join(wline) + '\n')

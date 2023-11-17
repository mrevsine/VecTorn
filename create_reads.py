import random
import os

random.seed(23)

fa_file_name = "chr22_14M_15M.fa"
fa_seq_name = "chr22:14000000-14999999"
fa_seq_len = 1000000
reads_1_name = "reads_1.fastq"
reads_2_name = "reads_2.fastq"

pos_snps = {}

for i in range(10):
    pos = int(random.random() * fa_seq_len)
    snps = {}
    mosaic_count = int(random.random() * 4) + 1
    freq_perc = [1] * mosaic_count
    total = 20 - mosaic_count
    for i in range(mosaic_count-1):
        to_give = int(random.random() * total)
        freq_perc[i] += to_give
        total -= to_give
    freq_perc[-1] += total
    bases = ['A','C','G','T']
    pos_counts = {}
    for i in range(mosaic_count):
        base = bases[int(random.random() * len(bases))]
        bases.remove(base)
        pos_counts[base] = freq_perc[i]
    pos_snps[pos] = pos_counts


temp_file_name = 'temp.fa'

os.system("rm -f " + temp_file_name)

positions = list(pos_snps.keys())
positions.sort()

for position in positions:
    os.system("samtools faidx " + fa_file_name + " " + fa_seq_name + ":" + str(position) + "-" + str(position) + " >> " + temp_file_name)

ref_bases = {}
ref_bases_file = open(temp_file_name)
for position in positions:
    ref_bases_file.readline()
    ref_bases[position] = ref_bases_file.readline().strip()

os.system("rm -f " + temp_file_name)


for i in range(20):
    vcf_iter_filename = 'temp_' + str(i) + '.vcf'
    vcf_iter_file = open(vcf_iter_filename, 'w')
    for position in positions:
        if 'A' in pos_snps[position]:
            b = 'A'
        elif 'C' in pos_snps[position]:
            b = 'C'
        elif 'G' in pos_snps[position]:
            b = 'G'
        else:
            b = 'T'
        if b != ref_bases[position]:
            to_write = fa_seq_name + '\t' + str(position) + '\t.\t' + ref_bases[position]
            to_write += ('\t' + b + '\t.\t.\t\n')
            vcf_iter_file.write(to_write)
        pos_snps[position][b] -= 1
        if pos_snps[position][b] == 0:
            pos_snps[position].pop(b)
    vcf_iter_file.close()
    gen_iter_pre = 'temp_gen_' + str(i)
    c = "perl simuG.pl -refseq " + fa_file_name + " -snp_vcf " + vcf_iter_filename
    c += (" -prefix " + gen_iter_pre)
    os.system(c)
    c = "./wgsim -1 150 -2 150 -r 0 -R 0 -X 0 -e 0.02 -N 3500 " + gen_iter_pre + ".simseq.genome.fa "
    c += (gen_iter_pre + "_1.fastq " + gen_iter_pre + "_2.fastq")
    os.system(c)
    os.system("cat " + gen_iter_pre + "_1.fastq >> " + reads_1_name)
    os.system("cat " + gen_iter_pre + "_2.fastq >> " + reads_2_name)
    os.system("rm -f " + gen_iter_pre + ".refseq2simseq.SNP.vcf")
    os.system("rm -f " + gen_iter_pre + ".refseq2simseq.map.txt")
    os.system("rm -f " + gen_iter_pre + ".simseq.genome.fa")
    os.system("rm -f " + gen_iter_pre + "_1.fastq")
    os.system("rm -f " + gen_iter_pre + "_2.fastq")
    os.system("rm -f " + vcf_iter_filename)



# Create [out].1.fastq and [out].2.fastq files containing 
# simulated paired-end reads from both a reference sequence
# and simulated recombinant sequences
# @authors Mahler Revsine, Sparsh Shah

###=============================================================================
# Imports

import math
import numpy as np
import os
import random


###=============================================================================
# Global parameters

# Seed params
sim_seed = 1
simuG_seed = 201812201903
wgsim_seed = 23

# Recombinant sequences simulation params
n_mutation_sims_per_recomb = 2

# Mosaic reference sequences simulation params
n_mosaic_sim_seqs = 20
n_mosaic_SNP_loci = 10

# Paired-end reads simulation params
wgsim_read_length = 150
wgsim_mutation_rate = 0
wgsim_indel_frac = 0
wgsim_indel_extension_prob = 0
wgsim_base_error_rate = 0.02

# Overall simulation params
ref_coverage = 250
recombinant_coverage = 1000

# Housekeeping params
temp_path = "data/mosaic"
temp_mosaic_sim_file_prefix = "temp_sim"
delete_intermediate_files = True


###=============================================================================
# Helper functions

# read in fasta file
# returns headers, sequences for the file
def read_fasta(fa_file):
  # read in file
  with open(fa_file, "r") as f:
    lines = f.readlines()
  lines = [line.rstrip() for line in lines]
  # parse lines, creating headers and sequences
  headers = []
  seqs = []
  curr_seq = ""
  for line in lines:
    if line.startswith(">"):
      if curr_seq != "":
        seqs.append(curr_seq)
        curr_seq = ""
      headers.append(line[1:])
    else:
      curr_seq += line
  if curr_seq != "":
    seqs.append(curr_seq)
  # return data
  return headers, seqs


# write sequences to fasta
def write_fasta(names, seqs, out_fa_file, max_char_per_line=80):
  with open(out_fa_file, "w") as f:
    for name, seq in zip(names, seqs):
      f.write(">" + name + "\n")
      for i in range(0, len(seq), max_char_per_line):
        f.write(seq[i:i+max_char_per_line] + "\n")


# Remove a polyA tail from an input sequence, if it exists
# Return cleaned sequence
# 4 A's would randomly occur ~ 0.25% of the time, seems like a good threshold 
def remove_polyA_tail(seq, min_tail_length=4):
  ending_A_count = 0
  for c in seq[::-1]:
    if c.upper() == "A":
      ending_A_count += 1
    else:
      break
  if ending_A_count < min_tail_length:
    return seq
  return seq[:-ending_A_count]


# Create a random SNP in a sequence
# Assume seed has been set previously
def create_random_SNP(seq):
  SNP_locus = random.randrange(len(seq))
  seq_locus_base = seq[SNP_locus].upper()
  SNP_base_options = [c for c in ["A","C","G","T"] if c != seq_locus_base]
  SNP_new_base = SNP_base_options[random.randrange(3)]
  SNP_seq = seq[:SNP_locus] + SNP_new_base + seq[SNP_locus+1:]
  return SNP_locus, SNP_new_base, SNP_seq
  

# Create a random SNP in each sequence from a list
def create_random_SNPs(seqs):
  SNP_loci = [random.randrange(len(seq)) for seq in seqs]
  SNP_bases = [random.randrange(3) for _ in range(len(seqs))]
  SNP_seqs = []
  for i in range(len(seqs)):
    seq = seqs[i]
    SNP_locus = SNP_loci[i]
    seq_locus_base = seq[SNP_locus].upper()
    SNP_base_options = [c for c in ["A","C","G","T"] if c != seq_locus_base]
    SNP_new_base = SNP_base_options[SNP_bases[i]]
    SNP_seqs.append(seq[:SNP_locus] + SNP_new_base + seq[SNP_locus+1:])
  return SNP_seqs


# choose locus within a vector Multiple Cloning Site at which to split
# for now just choose the midpoint
# in the future, we may randomize this selection process
def get_split_index(start_idx, end_idx):
  range = end_idx - start_idx
  half_range = int(range // 2) # integer division is analagous to floor()
  return start_idx + half_range


# create a vector insert sequence starting with vector seq, insert seq, and position at which to split the vector seq
def perform_insert(vector_seq, vector_split_pos, insert_seq):
  return vector_seq[:vector_split_pos] + insert_seq + vector_seq[vector_split_pos:]


# Starting with N genes and V vectors, create N*V recombinant sequences by inserting each gene into each vector
# For each recombinant seq, create a SNP at a random position within the inserted gene
def simulate_recombinant_sequences(gene_dict, vector_dict, vector_MCS_position_dict, out_fa):

  recomb_seqs = []
  recomb_seq_names = []
 
  # simulate for each pair of gene and vector
  for n,(gene_name,gene_seq) in enumerate(gene_dict.items()):
    for v,(vector_name,vector_seq) in enumerate(vector_dict.items()):

      # Simulate `n_mutation_sims_per_recomb` times per gene,vector pair
      # In Vecuum, this was 2 times
      last_locus = -1
      last_base = ""
      for i in range(n_mutation_sims_per_recomb):

        # Make sure that all simulated mutations per gene,vector pair are unique
        SNP_locus, SNP_base, SNP_gene = create_random_SNP(gene_seq)      
        while (SNP_locus == last_locus) and (SNP_base == last_base):
          SNP_locus, SNP_base, SNP_gene = create_random_SNP(gene_seq)      
        last_locus = SNP_locus
        last_base = SNP_base

        # Perform split
        vec_split_idx = get_split_index(*vector_MCS_position_dict[vector_name])
        recomb_seq = perform_insert(vector_seq, vec_split_idx, SNP_gene)
        recomb_seqs.append(recomb_seq)

        # Name is of the form vectorName:vectorSplitIdx_geneName:geneLength_SNPLocus:OriginalBase:MutatedBase
        recomb_seq_name = vector_name + ":" + str(vec_split_idx) + "_" + gene_name + ":" + str(len(SNP_gene)) + "_" + str(SNP_locus) + ":" + gene_seq[SNP_locus] + ":" + SNP_base
        recomb_seq_names.append(recomb_seq_name) 

  # Write recombinant sequences to file
  write_fasta(recomb_seq_names, recomb_seqs, out_fa)


# Create probability tables for SNPs at randomly selected loci
# returns dict of {position1:["A":p(A), "C":p(C), etc.], position2:[...}
def create_SNP_probabilities(n_loci, fa_seq_len):
  pos_snps = {} 
  max_prob = 20 # simuG treats each "point" of probability as 5%, so 20 = 100%
  for i in range(n_mosaic_SNP_loci):
    pos = int(random.random() * fa_seq_len)
    snps = {}
    mosaic_count = int(random.random() * 4) + 1 # rand in {1-4}
    freq_perc = [1] * mosaic_count
    total = max_prob - mosaic_count
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
  return pos_snps


# Get reference sequence base at a list of positions
# Returns dict of {position:base}
def get_reference_bases(ref_fa, ref_seq_name, positions):
  temp_file_name = temp_path + "/temp.fa"
  os.system("rm -f " + temp_file_name)
  os.system("touch " + temp_file_name)
  for position in positions:
    os.system("samtools faidx " + ref_fa + " " + ref_seq_name + ":" + str(position) + "-" + str(position) + " >> " + temp_file_name)
  ref_bases = {}
  ref_bases_file = open(temp_file_name)
  for position in positions:
    ref_bases_file.readline()
    ref_bases[position] = ref_bases_file.readline().strip()
  return ref_bases


# Create simulated mosaic reference sequences
def simulate_mosaic_ref_seqs(n_sim_seqs, n_loci, ref_fa, ref_seq_name, ref_seq_len, out_fa):
  
  # Create output directory
  os.system(f"rm -rf {temp_path}")
  os.makedirs(temp_path)

  # Get SNP probabilities
  pos_snps = create_SNP_probabilities(n_loci, ref_seq_len)

  # Get reference bases at SNP loci
  positions = list(pos_snps.keys())
  positions.sort()
  ref_bases = get_reference_bases(ref_fa, ref_seq_name, positions)

  # Simulate mosaic mutations of reference seq
  for i in range(n_sim_seqs):
    vcf_iter_filename = f"{temp_path}/{temp_mosaic_sim_file_prefix}_{i}.vcf"
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
        to_write = ref_seq_name + '\t' + str(position) + '\t.\t' + ref_bases[position]
        to_write += ('\t' + b + '\t.\t.\t\n')
        vcf_iter_file.write(to_write)
      pos_snps[position][b] -= 1
      if pos_snps[position][b] == 0:
        pos_snps[position].pop(b)
    vcf_iter_file.close()
    c = f"perl lib/simuG.pl -refseq {ref_fa} -snp_vcf {vcf_iter_filename} -seed {simuG_seed} -prefix {temp_path}/{temp_mosaic_sim_file_prefix}_{i}"
    os.system(c)

  # Concatenate all mosaic seqs into a single fasta
  os.system(f"touch {out_fa}")
  for i in range(n_sim_seqs):
    os.system(f"cat {temp_path}/{temp_mosaic_sim_file_prefix}_{i}.simseq.genome.fa >> {out_fa}")   

  # return SNP probability table in case user wants to save them
  return pos_snps


# Simulated taking paired-end reads from a fasta file using wgsim
def sim_paired_end_reads(fa_file, n_reads, out_temp_file_prefix):
  wgsim_cmd = f"wgsim -1 {wgsim_read_length} -2 {wgsim_read_length} -r {wgsim_mutation_rate} -R {wgsim_indel_frac} -X {wgsim_indel_extension_prob} -e {wgsim_base_error_rate} -N {n_reads} -S {wgsim_seed} {fa_file} {out_temp_file_prefix}.1.fastq {out_temp_file_prefix}.2.fastq"
  os.system(wgsim_cmd)


# Full pipeline
# Simulate recombinant reads and mosaic reference sequences
# Then take paired end reads from sequences
def simulate_reads_with_vector_contamination(out_path, gene_dict, vector_dict, vector_MCS_position_dict, ref_fa, ref_seq_name, ref_seq_length): 
  
  print("RUNNING FULL SIMULATION")

  # Create output directory
  if out_path.endswith("/"):
    out_path = out_path[:-1]
  if not os.path.isdir(out_path):
    os.makedirs(out_path)

  # Simulate recombinant sequences
  print("SIMULATING RECOMBINANT SEQUENCES")
  recomb_seqs_fa = out_path + "/recomb_seqs.fa"
  simulate_recombinant_sequences(gene_dict, vector_dict, vector_MCS_position_dict, recomb_seqs_fa)

  # Simulate mosaic sequences  
  print("SIMULATING MOSAIC SEQUENCES")
  mosaic_ref_seqs_fa = out_path + "/mosaic_ref_seqs.fa"
  simulate_mosaic_ref_seqs(n_mosaic_sim_seqs, n_mosaic_SNP_loci, ref_fa, ref_seq_name, ref_seq_length, mosaic_ref_seqs_fa)

  # Simulate paired-end reads from simulated recombinant sequences
  # First, get total length of recombinant sequences to determine num reads
  print("SIMULATING PAIRED-END READS FROM RECOMBINANT SEQUENCES")
  recomb_reads_prefix = temp_path + "/recomb_reads"
  total_recomb_seq_length = 0
  n_genes = len(gene_dict.keys())
  n_vectors = len(vector_dict.keys())
  for seq in gene_dict.values():
    total_recomb_seq_length += len(seq) * n_vectors
  for seq in vector_dict.values():
    total_recomb_seq_length += len(seq) * n_genes
  total_recomb_seq_length *= n_mutation_sims_per_recomb
  n_recomb_reads = int((total_recomb_seq_length * recombinant_coverage) / (2 * wgsim_read_length))
  print(f"Taking {n_recomb_reads} reads from {total_recomb_seq_length}-length recombinant genome")
  sim_paired_end_reads(recomb_seqs_fa, n_recomb_reads, recomb_reads_prefix)

  # Simulate paired-end reads from mosaic reference sequences
  # First, get total length of mosaic sequences to determine num reads
  print("SIMULATING PAIRED-END READS FROM MOSAIC REFERENCE SEQUENCES")
  mosaic_ref_reads_prefix = temp_path + "/mosaic_ref_reads"
  total_mosaic_seq_length = ref_seq_length * n_mosaic_sim_seqs
  n_mosaic_reads = int((total_mosaic_seq_length * ref_coverage) / (2 * wgsim_read_length))
  print(f"Taking {n_mosaic_reads} reads from {total_mosaic_seq_length}-length mosaic genome")
  sim_paired_end_reads(mosaic_ref_seqs_fa, n_mosaic_reads, mosaic_ref_reads_prefix)

  # Concatenate simulated reads into one file
  print("COMBINING ALL READS")
  out_fq1 = f"{out_path}/reads.1.fastq"
  out_fq2 = f"{out_path}/reads.2.fastq"
  os.system(f"touch {out_fq1}")
  os.system(f"touch {out_fq2}")
  os.system(f"cat {recomb_reads_prefix}.1.fastq >> {out_fq1}")
  os.system(f"cat {recomb_reads_prefix}.2.fastq >> {out_fq2}")
  os.system(f"cat {mosaic_ref_reads_prefix}.1.fastq >> {out_fq1}")
  os.system(f"cat {mosaic_ref_reads_prefix}.2.fastq >> {out_fq2}")

  # Optionally delete all intermediate files
  if delete_intermediate_files:
    print("DELETING TEMPORARY FILES")
    os.system(f"rm -rf {temp_path}")


###=============================================================================
# Read in needed data

def __read_needed_data__():

  # Needed paths
  data_path = "data/"
  gene_path = data_path + "genes/"
  vector_path = data_path + "vectors/"

  # Needed files
  gene_fa_files = [s for s in os.listdir(gene_path) if s.endswith(".fa")]
  vector_pos_bed = vector_path + "clone_vec.bed"

  # Read in gene fastas
  # Convert mRNA sequences to CDNA by removing polyA tail
  gene_CDNA_dict = {}
  for gene_fa_file in gene_fa_files:
    gene_name = gene_fa_file.split(".")[0]
    with open(gene_path + gene_fa_file, "r") as f:
      lines = f.readlines()
    mrna_seq = "".join([line.rstrip() for line in lines[1:]])
    cdna_seq = remove_polyA_tail(mrna_seq)
    gene_CDNA_dict[gene_name] = cdna_seq

  # read vector Multiple Cloning Site position file
  with open(vector_pos_bed, "r") as f:
    lines = f.readlines()
  vector_MCS_positions = {}
  for line in lines[1:]:
    line_terms = line.rstrip().split("\t")
    vector_MCS_positions[line_terms[0]] = (int(line_terms[1]), int(line_terms[2]))
  vectors = list(vector_MCS_positions.keys())

  # read vector sequence files
  vector_seqs = {}
  for vector in vectors:
    with open(vector_path + vector + ".fa") as f:
      lines = f.readlines()
    vector_seq = "".join([line.rstrip() for line in lines[1:]])
    vector_seqs[vector] = vector_seq

  # Return needed variables
  return gene_CDNA_dict, vector_seqs, vector_MCS_positions


###============================================================================
# Perform full simulation

def __main__():

  # Specify reference seq 
  reference_seq_fa_file = "data/ref/chr22_23M_24M.fa"
  reference_seq_name = "chr22:23000000-23999999"
  reference_seq_length = 1000000

  # Load needed data
  gene_dict, vector_dict, vector_MCS_positions_dict = __read_needed_data__()

  # For simplicity, use 10 genes and 10 vectors
  # Remove KP126802.1, the longest one
  del vector_MCS_positions_dict["KP126802.1"]
  del vector_dict["KP126802.1"]

  # Run full data simulation  
  out_path = "data/full_sim2"
  random.seed(sim_seed)
  simulate_reads_with_vector_contamination(out_path, gene_dict, vector_dict, vector_MCS_positions_dict, reference_seq_fa_file, reference_seq_name, reference_seq_length)


if __name__ == "__main__":
  __main__()




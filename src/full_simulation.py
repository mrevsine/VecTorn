# Code to produce simulated reads of vector inserts


###=============================================================================
# Imports

import math
import numpy as np
import os
import random


###=============================================================================
# Functions

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
def create_random_SNPs(seqs, seed=1):
  random.seed(seed)
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
def get_split_index(start_idx, end_idx, seed=1):
  #random.seed(seed) # not needed yet
  range = end_idx - start_idx
  half_range = int(range // 2) # integer division is analagous to floor()
  return start_idx + half_range


# create a vector insert sequence starting with vector seq, insert seq, and position at which to split the vector seq
def perform_insert(vector_seq, vector_split_pos, insert_seq):
  return vector_seq[:vector_split_pos] + insert_seq + vector_seq[vector_split_pos:]


# Starting with N genes and V vectors, create N*V recombinant sequences by inserting each gene into each vector
# For each recombinant seq, create a SNP at a random position within the inserted gene
def create_recombinant_sequences(gene_dict, vector_dict, vector_MCS_position_dict, n_mutation_sims_per_recomb=2, seed=1):

  random.seed(seed)  
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

  return recomb_seq_names, recomb_seqs


###=============================================================================# Read in needed data

# needed files
gene_dir = "data/genes/"
gene_fa_files = [s for s in os.listdir(gene_dir) if s.endswith(".fa")]
vector_dir = "data/vectors/"
vector_pos_bed = vector_dir + "clone_vec.bed"
reference_seq_fa_file = "data/ref/chr22_23M_24M.fa"

# Read in gene fastas
# Convert mRNA sequences to CDNA by removing polyA tail
gene_CDNA_dict = {}
for gene_fa_file in gene_fa_files:
  gene_name = gene_fa_file.split(".")[0]
  with open(gene_dir + gene_fa_file, "r") as f:
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
  with open(vector_dir + vector + ".fa") as f:
    lines = f.readlines()
  vector_seq = "".join([line.rstrip() for line in lines[1:]])
  vector_seqs[vector] = vector_seq


###============================================================================
# Create recombinant vectors

# For simplicity, use 10 genes and 10 vectors
# Remove KP126802.1, the longest one
del vector_MCS_positions["KP126802.1"]
del vector_seqs["KP126802.1"]
vectors = [v for v in vectors if v != "KP126802.1"]

# Create recombinant sequences with mutations
recomb_names, recomb_seqs = create_recombinant_sequences(gene_CDNA_dict, vector_seqs, vector_MCS_positions)

# Write recombinant seqs to file
out_fa_file = "data/simulated_reads/recomb_vectors_b.fa"
write_fasta(recomb_names, recomb_seqs, out_fa_file)




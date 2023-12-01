# Code to produce simulated reads of vector inserts

import math
import numpy as np
import random

###=============================================================================# Functions

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


# sample sequence(s) from the supplied full_seq sequence
# samples are generated under a Gaussian model of sample length with the mean_length and sd_length parameters. If no sd_length is given, it will be set to the sqrt of the mean_length
# then for each length sample, its position is randomly selected under a uniform model across the full sequence
# defaults to returning a single sample, but can return any number specified by n_samples
def get_insert_seq_params(full_seq_length, mean_length=10000, sd_length=None, n_samples=1, seed=1):
  
  # process parameters
  if n_samples < 1:
    print("Invalid number of requested samples. Must be at least 1.")
    exit()
  if sd_length is None:
    sd_length = math.sqrt(mean_length) # poisson assumption

  # set seed
  np.random.seed(seed)

  # randomly generate lengths of each sample
  # any lengths longer than the input sequence will be overwritten as the length of the input sequence
  sample_lengths = np.random.normal(loc = mean_length, scale = sd_length, size = n_samples)
  sample_lengths = [int(l) for l in sample_lengths]
  sample_lengths = [max(0,l) for l in sample_lengths]
  sample_lengths = [min(l, full_seq_length) for l in sample_lengths]
  
  # randomly generate positions of each sample
  sample_positions = np.array([np.random.randint(0, full_seq_length-l+1) for l in sample_lengths])

  # return [(position_i, length_i)] tuple array for sampled sequences
  return [(p,l) for p,l in zip(sample_positions, sample_lengths)]


# select a vector among options
def choose_vectors(vector_options, n_samples=1, seed=1):
  random.seed(seed)
  return random.choices(vector_options, k=n_samples)


# choose locus within a vector Multiple Cloning Site at which to split
# for now just choose the midpoint
# in the future, we may randomize this selection process
def get_split_index(start_idx, end_idx, seed=1):
  # seed(seed) # not needed yet
  range = end_idx - start_idx
  half_range = int(range // 2) # integer division is analagous to floor()
  return start_idx + half_range


# select vector insertion details n_samples times
# return list of tuples of (vector name, vector insertion position) for each sample
def get_vector_insertion_params(vector_options, vector_positions, n_samples=1, all_same_vector=False, seed=1):
  
  # choose vectors and their split positions
  if all_same_vector:
    vector_choice = choose_vectors(vector_options, n_samples=1, seed=seed)[0]
    vector_choices = [vector_choice for _ in range(n_samples)]
  else:
    vector_choices = choose_vectors(vector_options, n_samples=n_samples, seed=seed)

  # for each vector, get insertion position
  vector_inserts = [get_split_index(*vector_positions[vec], seed=seed) for vec in vector_choices]

  # return tuple list of [(vector_i, split_index_i)]
  return [(c,i) for c,i in zip(vector_choices, vector_inserts)]


# create a vector insert sequence starting with vector seq, insert seq, and position at which to split the vector seq
def perform_insert(vector_seq, vector_split_pos, insert_seq):
  return vector_seq[:vector_split_pos] + insert_seq + vector_seq[vector_split_pos:]


# create a vector insert sequence, starting with no needed parameters
# uses rng to select several parameters
def get_insertions(full_seq, vector_names, vector_seqs, vector_positions, n_samples=1, mean_insert_length=10000, sd_insert_length=None, all_same_vector=False, ref_name="ref", seed=1):

  # generate sample insertion sequence params from full sequence
  insert_seq_params = get_insert_seq_params(len(full_seq), \
                                            mean_length=mean_insert_length, \
                                            sd_length = sd_insert_length, \
                                            n_samples=n_samples, \
                                            seed=seed)

  # generate vector insert data for each insertion sequence
  vector_insert_params = get_vector_insertion_params(vector_options=vector_names, vector_positions=vector_positions, n_samples=n_samples, all_same_vector=all_same_vector, seed=seed)  

  # get recombinant vector insertion sequences
  recombinant_seqs = [perform_insert(vector_seqs[vec_name], split_idx, full_seq[insert_idx:insert_idx+insert_length]) for ((vec_name, split_idx), (insert_idx, insert_length)) in zip(vector_insert_params, insert_seq_params)]

  # get recombinant insert names for the fasta headers
  recombinant_names = [ref_name + ":" + str(insert_idx) + "-" + str(insert_idx+insert_length) + "_" + vec_name + ":" + str(split_idx) for ((vec_name, split_idx), (insert_idx, insert_length)) in zip(vector_insert_params, insert_seq_params)]

  # return recombinant names and sequences
  return recombinant_names, recombinant_seqs


###=============================================================================# Read in needed data

# needed files
vector_dir = "data/vectors/"
vector_pos_bed = vector_dir + "clone_vec.bed"
reference_seq_fa_file = "data/ref/chr22_14M_15M.fa"

# read vector position file
with open(vector_pos_bed, "r") as f:
  lines = f.readlines()
vector_positions = {}
for line in lines[1:]:
  line_terms = line.rstrip().split("\t")
  vector_positions[line_terms[0]] = (int(line_terms[1]), int(line_terms[2]))
vectors = list(vector_positions.keys())

# read vector sequence files
vector_seqs = {}
for vector in vectors:
  with open(vector_dir + vector + ".fa") as f:
    lines = f.readlines()
  vector_seq = "".join([line.rstrip() for line in lines[1:]])
  vector_seqs[vector] = vector_seq


###=============================================================================# Test functionality

run_tests = False

if run_tests:

  ## read_fasta(fa_file)
  print("Testing read_fasta()")
  ref_headers, ref_seqs = read_fasta(reference_seq_fa_file)
  print("%s headers, %s seqs" % (len(ref_headers), len(ref_seqs)))
  ref_header = ref_headers[0]
  ref_seq = ref_seqs[0]
  print(ref_header)
  print(len(ref_seq))
  print("--------------------")

  ## get_insert_seq_params(full_seq_length, mean_length=10000, sd_length=None, n_samples=1, seed=1):
  print("Testing get_insert_seq_params()")
  insert_seq_params = get_insert_seq_params(len(ref_seq), n_samples = 10)
  print(insert_seq_params)
  insert_seq_params = get_insert_seq_params(len(ref_seq), n_samples = 10)
  print(insert_seq_params)
  print("--------------------")

  ## choose_vectors(vector_options, n_samples=1, seed=1):
  print("Testing choose_vectors()")
  vector_choices = choose_vectors(vectors, n_samples = 10)
  print(vector_choices)
  vector_choices = choose_vectors(vectors, n_samples = 10)
  print(vector_choices)
  print("--------------------")

  ## get_split_index(start_idx, end_idx, seed=1):
  print("Testing get_split_index()")
  for k,v in vector_positions.items():
    print("%s: split of %s and %s = %s" % (k, v[0], v[1], get_split_index(v[0],v[1])))
  print("--------------------")

  ## get_vector_insertion_params(vector_options, vector_positions, n_samples=1, all_same_vector=False, seed=1):
  print("Testing get_vector_insertion_params()")
  vector_insertion_params = get_vector_insertion_params(vectors, vector_positions, n_samples = 10, all_same_vector = False)
  print(vector_insertion_params)
  vector_insertion_params = get_vector_insertion_params(vectors, vector_positions, n_samples = 10, all_same_vector = False)
  print(vector_insertion_params)
  vector_insertion_params = get_vector_insertion_params(vectors, vector_positions, n_samples = 10, all_same_vector = True)
  print(vector_insertion_params)
  print("--------------------")

  ## get_insertions(full_seq, ref_name="ref", vector_names, vector_seqs, vector_positions, n_samples=1, mean_insert_length=10000, sd_insert_length=None, all_same_vector=False, seed=1):
  print("Testing perform_insert(), get_insertions(), and write_fasta()")
  recomb_names, recomb_seqs = get_insertions(full_seq = ref_seq, \
                                             vector_names = vectors,
                                             vector_seqs = vector_seqs,
                                             vector_positions =vector_positions, \
                                             n_samples = 10,
                                             seed = 1)
  print("Writing fasta")
  results_fa_file = "results/test_results_1.fa"
  write_fasta(recomb_names, recomb_seqs, results_fa_file)

  # Now read it back in and check that all data matches
  print("Testing fasta")
  in_recomb_names, in_recomb_seqs = read_fasta(results_fa_file)
  print("Do headers match?", recomb_names == in_recomb_names)
  print("Do sequences match?", recomb_seqs == in_recomb_seqs)
  for header,seq in zip(in_recomb_names, in_recomb_seqs):
    ref_start = int(header.split("_")[0].split(":")[1].split("-")[0])
    ref_end = int(header.split("_")[0].split(":")[1].split("-")[1])
    ref_seq_segment = ref_seq[ref_start:ref_end]
    vector_name = header.split("_")[1].split(":")[0]
    vector_split = int(header.split("_")[1].split(":")[1])
    vector_seq_1 = vector_seqs[vector_name][:vector_split]
    vector_seq_2 = vector_seqs[vector_name][vector_split:]
    full_check_seq = vector_seq_1 + ref_seq_segment + vector_seq_2
    print("Does " + header + " seq match?", full_check_seq == seq)
    if full_check_seq != seq:
      print("Expected length is %s but observed length is %s" % (len(seq), len(full_check_seq)))
  print("--------------------")

  print("All tests done")

###=============================================================================
# Use script

ref_headers, ref_seqs = read_fasta(reference_seq_fa_file)
ref_seq = ref_seqs[0]
run_seed = 123
run_n_samples = 100
run_insert_length = 10000
run_out_fa_file = "data/simulated_reads/100_samples_10k_length_123_seed_b.fa"
#run_out_fa_file = "results/" + str(run_n_samples) + "samples_" + str(run_seed) + "seed.fa"
recomb_names, recomb_seqs = get_insertions(full_seq = ref_seq, \
                                           vector_names = vectors,
                                           vector_seqs = vector_seqs,
                                           vector_positions =vector_positions, \
                                           n_samples = run_n_samples,
                                           mean_insert_length = run_insert_length, \
                                           seed = run_seed)
write_fasta(recomb_names, recomb_seqs, run_out_fa_file)


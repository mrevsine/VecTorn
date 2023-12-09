# Detect vector contamination using a HMM
# For now, model inputs are a kmerized sequence
# Learn kmer emissions in human and vector
# Use this to predict sites of vector contamination

###=============================================================================
# Imports

import itertools
from matplotlib import pyplot as plt
import numpy as np
import os
import random

###=============================================================================
# Set seed

np_seed = 42
np.random.seed(np_seed)

###=============================================================================
# My implementation of an HMM

class MyHMM:

  # emissions is an n_component x n_feature np matrix
  # transitions is an n_component x n_component np matrix
  # starting states is an n_component x 1 np matrix
  def __init__(self, n_components, emissions, transitions, starting_states):
    self.n_components = n_components
    self.E = np.log(emissions)
    self.T = np.log(transitions)
    self.S = np.log(starting_states)

  # X is a single sample
  def viterbi(self,X):

    # Initialize needed tables
    prediction_table = np.zeros((self.n_components, len(X)))
    backtrack_table = np.zeros((self.n_components, len(X)), dtype=int)
    
    # Full in the first column
    prediction_table[:,0] = self.S + self.E[:, X[0]]
    backtrack_table[:,0] = [i for i in range(self.n_components)]

    # Compute the rest of the columns
    for i in range(1, len(X)):
      for s in range(self.n_components):
        transition_probs = prediction_table[:,i-1] + self.T[s,:]
        max_trans_prob = np.max(transition_probs)
        prediction_table[s,i] = max_trans_prob + self.E[s,X[i]]
        backtrack_table[s,i] = np.argmax(transition_probs)

    # Get optimal path
    optimal_path = [np.argmax(prediction_table[:,-1])]
    for i in range(len(X)-1, 0, -1):
      optimal_path.append(backtrack_table[optimal_path[-1], i])
    optimal_path.reverse()

    # Return best path and its corresponding log probability
    return optimal_path, np.max(prediction_table[:,-1])

  # X is list of samples
  # Returns a list of (optimal_path, log prob) tuples
  def predict(self,X):
    if isinstance(X[0], list): # if multiple samples
      predictions = [self.viterbi(X[i]) for i in range(len(X))]
    else:
      predictions = [self.viterbi(X)]
    return predictions

###=============================================================================
# Functions

# Return list of all possbile kmers of length `k`
def get_all_possible_kmers(k):
  bases = ["A","C","G","T"]
  kmers = [''.join(seq) for seq in itertools.product(bases, repeat=k)]
  return kmers

# Convert DNA sequence into encoded kmer sequence
# Encoded by representing each kmer as an integer specifying its alphabetical position among all possible kmers
def kmerize(seq, k, overlapping=False):
  possible_kmers = get_all_possible_kmers(k)
  kmer_encoding_table = {seq:i for i,seq in enumerate(possible_kmers)}
  if overlapping:
    encoding = [kmer_encoding_table[seq[i:i+k]] for i in range(0,len(seq)-k+1)]
  else:
    end_index = len(seq) - (len(seq) % k) # cap seq length at a multiple of k 
    encoding = [kmer_encoding_table[seq[i:i+k]] for i in range(0,end_index, k)]
  return encoding

# Get kmer frequencies within a sequence
# Returns a list ordered by alphabetical position of each possible kmer
# Each element is the fraction of kmers within the sequence matching that kmer
def get_kmer_frequencies(seq, k, overlapping=False, revcomp=True):
  seq_kmers = kmerize(seq, k, overlapping)
  if revcomp:
    seq_kmers += kmerize(get_revcomp(seq), k, overlapping)
  n_seq_kmers = float(len(seq_kmers))
  kmer_counts = [0 for _ in range(4 ** k)]
  for kmer in seq_kmers:
    kmer_counts[kmer] += 1
  kmer_frequencies = [ct / n_seq_kmers for ct in kmer_counts]
  return kmer_frequencies

# Given a list of DNA sequences, builds an emission probability table
# Each sequence is considered to represent a unique component
# The table is len(seqs) x 4^k
def build_emission_probability_table(seqs, k, overlapping=False, revcomp=True):
  E = np.zeros((len(seqs), 4**k))
  for i,seq in enumerate(seqs):
    seq_kmer_frequencies = get_kmer_frequencies(seq, k, overlapping, revcomp)
    E[i,:] = seq_kmer_frequencies
  return E

# Return a "clean" version of the sequence where every base is in {A,C,G,T}
# Convert lowercase to uppercase and discard any wildcard characters
def clean_seq(seq):
  seq = seq.upper()
  seq = "".join([c for c in seq if c in ["A","C","G","T"]])
  return seq

# Reverse complement of a sequence
def get_revcomp(seq):
  revcomp_table = {"A":"T","C":"G","G":"C","T":"A"}
  revcomp_seq = "".join([revcomp_table[c] for c in seq])[::-1]
  return revcomp_seq

# Return visualization of emission probability table
# Bar chart of probabilities for each component
# Compare any 2 rows of emission table
def visualize_emission_probability_table(E, i1=0, i2=1, name1="Human", name2="Vector", relative=False, savename=None):
  if relative:
    fig = plt.figure()
    plt.bar([i for i in range(E.shape[1])], E[i1,:]-E[i2,:])
    plt.ylabel("Relative emission probability")
    plt.title("Human - vector emission probability")
  else:
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.bar([i for i in range(E.shape[1])], E[i1,:])
    ax1.set_title(name1)
    ax1.set_ylabel("Emission probability")
    ax2.bar([i for i in range(E.shape[1])], E[i2,:])
    ax2.set_title(name2)
    ax2.set_ylabel("Emission probability")
  if savename:
    plt.savefig(savename)
  else:
    plt.show()
 
# Naively builds a transition probability table between states
# Assigns change of transitioning between any states to transition prob
def build_transition_probability_table(transition_prob, n_components):
  stay_prob = 1. - ((n_components-1) * transition_prob)
  T = np.eye(n_components, n_components) * stay_prob + (1 - np.eye(n_components, n_components)) * transition_prob
  return T 

# Initialize a CategoricalHMM from the training sequences
def get_HMM(training_seqs, k, overlapping=False, transition_prob=0.01, revcomp=True, visualize_emissions=False, savename=None):
  print("BUILDING MODEL WITH k=" + str(k))
  S = np.array([1 / float(len(training_seqs)) for _ in range(len(training_seqs))])
  T = build_transition_probability_table(transition_prob, len(training_seqs))
  E = build_emission_probability_table(training_seqs, k, overlapping, revcomp)
  
  if visualize_emissions:
    visualize_emission_probability_table(E, savename=savename)

  model = MyHMM(n_components=len(training_seqs), emissions=E, transitions=T, starting_states=S)
  return model

# Return a confusion matrix for MyHMM model output
# `label_list` is a list of labels
# `predictions` is a list of (backtrack, prob) values output by MyHMM.predict()
# For now, assume 2 classes
# Binary: whether to test all kmers (False) or treat a prediction as a vector based on whether any sequence is vector (True)
def confusion_matrix(label_list, prediction_results, n_classes=2, binary=True):
  conf_mat = np.zeros((n_classes, n_classes))
  for i in range(len(prediction_results)):
    predictions = prediction_results[i][0]
    labels = label_list[i]
    if binary:
      # if any prediction is vector, the sequence is vector
      yh = 0 if 0 in predictions else 1 
      y = 0 if 0 in labels else 1
      conf_mat[yh,y] += 1
    else:
      for yh,y in zip(predictions, labels):
        conf_mat[yh,y] += 1
  return conf_mat

# Get precision, recall, and F1 from confusion matrix
# For now, assume confusion matrix is 2x2
def conf_mat_stats(conf_mat):
  TP = conf_mat[0,0]
  FN = conf_mat[1,0]
  FP = conf_mat[0,1]
  TN = conf_mat[1,1]
  precision = TP / float(TP + FP)
  recall = TP / float(TP + FN)
  F1 = (2 * TP) / float((2 * TP) + FP + FN)
  return precision, recall, F1
   
###=============================================================================
# Load data

# Needed files
# We need 2 pieces of data - the reference human sequence and vector sequences
data_path = "data/"
ref_human_fa = data_path + "ref/chr22_23M_24M.fa"
ref_vector_fa = data_path + "univec/univec.fa"
simulation_path = data_path + "full_sim1/"
reads1_fq = simulation_path + "reads.1.fastq"
reads2_fq = simulation_path + "reads.2.fastq"

# Read in human region
with open(ref_human_fa, "r") as f:
  lines = f.readlines()
lines = [line.rstrip() for line in lines]
ref_human_seq = "".join([line for line in lines if not line.startswith(">")])
ref_human_seq = clean_seq(ref_human_seq)

# Read in vector sequences
# Concatenate them all into 1 sequence for now
with open(ref_vector_fa, "r") as f:
  lines = f.readlines()
lines = [line.rstrip() for line in lines]
ref_vector_seq = "".join([line for line in lines if not line.startswith(">")])
ref_vector_seq = clean_seq(ref_vector_seq)
del lines

"""
# Visualize emission probabilities
out_path = "results/HMM/figs/revcomp/"
training_seqs = [ref_human_seq, ref_vector_seq]
E_overlapping = build_emission_probability_table(training_seqs, k=5, overlapping=True)
E_non_overlapping = build_emission_probability_table(training_seqs, k=5, overlapping=False)
ot = visualize_emission_probability_table(E_overlapping, relative = True, savename = out_path + "overlapping_relative_5mer_frequency.png")
of = visualize_emission_probability_table(E_overlapping, relative = False, savename = out_path + "overlapping_5mer_frequency.png")
nt = visualize_emission_probability_table(E_non_overlapping, relative = True, savename = out_path + "non_overlapping_relative_5mer_frequency.png")
nf = visualize_emission_probability_table(E_non_overlapping, relative = False, savename = out_path + "non_overlapping_5mer_frequency.png")
exit()
"""

# Read in reads files
sim_headers_1 = []
sim_reads_1 = []
with open(reads1_fq, "r") as f:
  for i,line in enumerate(f):
    if i % 4 == 0:
      sim_headers_1.append(line.rstrip())
    elif i % 4 == 1:
      sim_reads_1.append(clean_seq(line.rstrip()))
sim_headers_2 = []
sim_reads_2 = []
with open(reads2_fq, "r") as f:
  for i,line in enumerate(f):
    if i % 4 == 0:
      sim_headers_2.append(line.rstrip())
    elif i % 4 == 1:
      sim_reads_2.append(clean_seq(line.rstrip()))

# Optionally sample N lines
n_sample_reads = 10000
random.seed(1)
sample_line_idxs = random.sample(range(0,len(sim_headers_1)), n_sample_reads)
sim_headers_1 = [sim_headers_1[i] for i in sample_line_idxs]
sim_headers_2 = [sim_headers_2[i] for i in sample_line_idxs]
sim_reads_1 = [sim_reads_1[i] for i in sample_line_idxs]
sim_reads_2 = [sim_reads_2[i] for i in sample_line_idxs]

# Label sequences
# Human is 1 (negative), vector is 0 (positive)
sim_seq_labels = [int(header.startswith("@chr22")) for header in sim_headers_1]
del sim_headers_1
del sim_headers_2

print("Loaded", len(sim_seq_labels), "simulated sequences")
print(str(len(sim_seq_labels)-sum(sim_seq_labels)) + " vector, " + str(sum(sim_seq_labels)) + " human") 

###=============================================================================
# Set revcomp for all trials

revcomp = True

###=============================================================================
# Test k vs. transition probability
# Hold overlap_train, overlap_test, and revcomp constant 

print("##### TESTING k vs. TRANSITION PROBABILITY ON SIMULATED DATA SEQUENCES #####")

training_seqs = [ref_vector_seq, ref_human_seq]
overlapping_train = False
overlapping_test = False

print("Global parameters:")
print("revcomp: " + str(revcomp) + ", overlapping_train: " + str(overlapping_train) + ", overlapping_test: " + str(overlapping_test))

t_probs = [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1]
ks = [k for k in range(1,9)]

S = np.array([1 / float(len(training_seqs)) for _ in range(len(training_seqs))])

for k in ks:

  E = build_emission_probability_table(training_seqs, k, overlapping_train, revcomp=revcomp)  
  testing_kmers = [kmerize(seq1,k,overlapping_test) + kmerize(seq2,k,overlapping_test) for seq1,seq2 in zip(sim_reads_1,sim_reads_2)]
  testing_labels = [[sim_seq_labels[i] for _ in range(len(testing_kmers[i]))] for i in range(len(testing_kmers))]

  for t_prob in t_probs:

    print("TRAINING ON TRANSITION PROBABILITY =", str(t_prob), "k =", str(k))
    T = build_transition_probability_table(t_prob, len(training_seqs))
    model = MyHMM(n_components=len(training_seqs), emissions=E, transitions=T, starting_states=S)

    predictions = model.predict(testing_kmers)

    conf_mat = confusion_matrix(testing_labels, predictions, binary=True)
    print("       ",conf_mat)
    stats = conf_mat_stats(conf_mat)
    print("        Precision = %s" % stats[0])
    print("        Recall = %s" % stats[1])
    print("        F1 = %s" % stats[2])

print("==============================================================")

###=============================================================================
# Test out model architectures on training dataset

print("##### TESTING k AND OVERLAP ON TRAINING DATA KMERS #####")

# Put sequences into a list
# 0th element is vector, 1st is human
training_seqs = [ref_vector_seq, ref_human_seq]

# Test on a variety of parameter values
ks = [k for k in range(1,9)]
overlapping_train_vals = [True,False]
overlapping_test_vals = [True,False]
for k in ks:
  for overlapping_train in overlapping_train_vals:
    print("TRAINING WITH k =", str(k), "AND OVERLAPPING =", overlapping_train)
    model = get_HMM(training_seqs, k, overlapping=overlapping_train, revcomp=revcomp)    
    for overlapping_test in overlapping_test_vals:
      print("    TESTING WITH OVERLAPPING =", overlapping_test)
      testing_kmers = [kmerize(seq,k,overlapping_test) for seq in training_seqs]
      testing_labels = [[i for _ in testing_kmers[i]] for i in range(len(testing_kmers))]
      predictions = model.predict(testing_kmers)
      conf_mat = confusion_matrix(testing_labels, predictions, binary=False)
      print("       ", conf_mat)
      stats = conf_mat_stats(conf_mat)
      print("        Precision = %s" % stats[0])
      print("        Recall = %s" % stats[1])
      print("        F1 = %s" % stats[2])

print("==============================================================")

###=============================================================================
# Build model

# Put sequences into a list
training_seqs = [ref_vector_seq, ref_human_seq]

# Initialize model
k = 3
overlapping_train = False
model = get_HMM(training_seqs, k, overlapping=overlapping_train, revcomp=revcomp)

overlapping_test = True
testing_kmers = [kmerize(seq,k,overlapping_test) for seq in training_seqs]
testing_labels = [[i for _ in testing_kmers[i]] for i in range(len(testing_kmers))]
predictions = model.predict(testing_kmers)
conf_mat = confusion_matrix(testing_labels, predictions)
print(conf_mat)
stats = conf_mat_stats(conf_mat)
print("Precision = %s" % stats[0])
print("Recall = %s" % stats[1])
print("F1 = %s" % stats[2])

###=============================================================================
# Test model on simulated dataset

# Get testing data
overlapping_test = False
overlap_str = "_o" if overlapping_test else "_u"
test_data_file = "results/" + str(k) + "mers" + overlap_str + ".txt"
if os.path.isfile(test_data_file):
  #print("LOADING KMERIZED TEST DATA")
  with open(test_data_file, "r") as f:
    lines = f.readlines()
  testing_kmers = [line.rstrip().split("\t") for line in lines]
  del lines
  testing_kmers = [[int(s) for s in terms] for terms in testing_kmers]
else:
  #print("KMERIZING TEST DATA")
  testing_kmers = [kmerize(seq1,k,overlapping_test) + kmerize(seq2,k,overlapping_test) for seq1,seq2 in zip(sim_reads_1,sim_reads_2)]
  print("SAVING TEST DATA")
  with open(test_data_file, "w") as f:
    for X in testing_kmers:
      f.write("\t".join([str(n) for n in X]) + "\n")
del sim_reads_1
del sim_reads_2

# For now, label all kmers in a sequence with their overall label
# i.e. if a read has any vector, call it all vector
testing_labels = [[sim_seq_labels[i] for _ in range(len(testing_kmers[i]))] for i in range(len(testing_kmers))]

# Test model
predictions = model.predict(testing_kmers)

# Get model accuracy
conf_mat = confusion_matrix(testing_labels, predictions)
print(conf_mat)
stats = conf_mat_stats(conf_mat)
print("Precision = %s" % stats[0])
print("Recall = %s" % stats[1])
print("F1 = %s" % stats[2])

###=============================================================================
# Test out model architectures on simulated reads

print("##### TESTING k AND OVERLAP ON SIMULATED DATA KMERS #####")

# Put sequences into a list
training_seqs = [ref_vector_seq, ref_human_seq]

# Test on a variety of parameter values
ks = [k for k in range(1,9)]
overlapping_train_vals = [True,False]
overlapping_test_vals = [True,False]
for k in ks:
  for overlapping_train in overlapping_train_vals:
    print("TRAINING WITH k =", str(k), "AND OVERLAPPING =", overlapping_train)
    model = get_HMM(training_seqs, k, overlapping=overlapping_train, revcomp=revcomp)    
    for overlapping_test in overlapping_test_vals:
      print("    TESTING WITH OVERLAPPING =", overlapping_test)

      testing_kmers = [kmerize(seq1,k,overlapping_test) + kmerize(seq2,k,overlapping_test) for seq1,seq2 in zip(sim_reads_1,sim_reads_2)]
      testing_labels = [[sim_seq_labels[i] for _ in range(len(testing_kmers[i]))] for i in range(len(testing_kmers))]

      predictions = model.predict(testing_kmers)

      conf_mat = confusion_matrix(testing_labels, predictions, binary=False)
      print("       ",conf_mat)
      stats = conf_mat_stats(conf_mat)
      print("        Precision = %s" % stats[0])
      print("        Recall = %s" % stats[1])
      print("        F1 = %s" % stats[2])

print("==============================================================")

###=============================================================================
# Test out model architectures on simulated reads

print("##### TESTING k AND OVERLAP ON SIMULATED DATA SEQUENCES #####")

# Put sequences into a list
training_seqs = [ref_vector_seq, ref_human_seq]

# Test on a variety of parameter values
ks = [k for k in range(1,9)]
overlapping_train_vals = [True,False]
overlapping_test_vals = [True,False]
for k in ks:
  for overlapping_train in overlapping_train_vals:
    print("TRAINING WITH k =", str(k), "AND OVERLAPPING =", overlapping_train)
    model = get_HMM(training_seqs, k, overlapping=overlapping_train, revcomp=revcomp)    
    for overlapping_test in overlapping_test_vals:
      print("    TESTING WITH OVERLAPPING =", overlapping_test)

      testing_kmers = [kmerize(seq1,k,overlapping_test) + kmerize(seq2,k,overlapping_test) for seq1,seq2 in zip(sim_reads_1,sim_reads_2)]
      testing_labels = [[sim_seq_labels[i] for _ in range(len(testing_kmers[i]))] for i in range(len(testing_kmers))]

      predictions = model.predict(testing_kmers)

      conf_mat = confusion_matrix(testing_labels, predictions, binary=True)
      print("       ",conf_mat)
      stats = conf_mat_stats(conf_mat)
      print("        Precision = %s" % stats[0])
      print("        Recall = %s" % stats[1])
      print("        F1 = %s" % stats[2])

print("==============================================================")

###=============================================================================
# For one successful model architecture, test effect of transition_prob

print("##### TESTING TRANSITION PROBABILITY ON SIMULATED DATA SEQUENCES #####")

training_seqs = [ref_vector_seq, ref_human_seq]
t_probs = [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1]
k = 5
overlapping_train = False
S = np.array([1 / float(len(training_seqs)) for _ in range(len(training_seqs))])
E = build_emission_probability_table(training_seqs, k, overlapping_train, revcomp=revcomp)
  
for t_prob in t_probs:
  print("TRAINING ON TRANSITION PROBABILITY =", str(t_prob), "k=", str(k), "OVERLAPPING =", overlapping_train)
  T = build_transition_probability_table(t_prob, len(training_seqs))
  model = MyHMM(n_components=len(training_seqs), emissions=E, transitions=T, starting_states=S)

  for overlapping_test in [True,False]:
    print("    TESTING ON OVERLAPPING =", overlapping_test)
    
    testing_kmers = [kmerize(seq1,k,overlapping_test) + kmerize(seq2,k,overlapping_test) for seq1,seq2 in zip(sim_reads_1,sim_reads_2)]
    testing_labels = [[sim_seq_labels[i] for _ in range(len(testing_kmers[i]))] for i in range(len(testing_kmers))]

    predictions = model.predict(testing_kmers)

    conf_mat = confusion_matrix(testing_labels, predictions, binary=True)
    print("       ",conf_mat)
    stats = conf_mat_stats(conf_mat)
    print("        Precision = %s" % stats[0])
    print("        Recall = %s" % stats[1])
    print("        F1 = %s" % stats[2])

print("==============================================================")

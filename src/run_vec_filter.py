import os
import subprocess
import numpy as np
import re

def calculate_results(true_labels, predicted_labels_1, predicted_labels_2):
    tp = np.sum((true_labels == 1) & ((predicted_labels_1 == 1) | (predicted_labels_2 == 1)))
    tn = np.sum((true_labels == 0) & (predicted_labels_1 == 0) & (predicted_labels_2 == 0))
    fp = np.sum((true_labels == 0) & ((predicted_labels_1 == 1) | (predicted_labels_2 == 1)))
    fn = np.sum((true_labels == 1) & (predicted_labels_1 == 0) & (predicted_labels_2 == 0))

    return tp, tn, fp, fn

def calculate_metrics(tp, tn, fp, fn):
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * (precision * recall) / (precision + recall)
    accuracy = (tp + tn) / len(true_labels)

    return precision, recall, f1, accuracy

def parse_log(log_file):
	time = 0.0
	mem = 0
	with open(log_file, 'r') as log_fp:
		for line in log_fp:
			time_match = re.match(r'Wall Time:\s+([\d.]+)', line)
			memory_match = re.match(r'Max Memory:\s+(\d+)', line)

			# Process matches
			if time_match:
				time = float(time_match.group(1))
			elif memory_match:
				mem = int(memory_match.group(1))

	return time, mem

def run_cmd (cmd):
	print(cmd)
	os.system(cmd)

def get_output_cmd (cmd):
	print(cmd)
	return subprocess.getoutput(cmd)

# ref_name = "data/ref/chr22"
vec_name = "data/univec/univec"

bf_run_src = "src/bf/bf_run.cpp"
bf_build_path = "build"
bf_run_prg = f"{bf_build_path}/bf_run"

reads_1_file = "data/full_sim/reads.1.fastq"
bf_result_1_file = "results/bf/bf_class_1.txt"
reads_2_file = "data/full_sim/reads.2.fastq"
bf_result_2_file = "results/bf/bf_class_2.txt"
reads_class = "data/full_sim/reads.class"

all_results_file = "results/bf/full_sim.txt"
log_file = "results/log"

if not os.path.exists(bf_build_path):
	os.makedirs(f"{bf_build_path}/")
if not os.path.isfile(bf_run_prg):
	os.system(f"g++ {bf_run_src} -o {bf_run_prg}")

if not os.path.isfile(f"{reads_class}"):
	run_cmd("awk '/^@/ {print /^@chr22/ ? 0 : 1}' " + f"{reads_1_file} > {reads_class}")

time_cmd = "/usr/bin/time --format=\"Wall Time: %e\nMax Memory: %M\""

k = 21
hit = 4

true_labels = np.loadtxt(reads_class, dtype=int)

run_cmd(f"{time_cmd} {bf_run_prg} {reads_1_file} {vec_name}.bf {k} {hit} {bf_result_1_file} 2> {log_file}")
time_1, mem_1 = parse_log(log_file)
run_cmd(f"{time_cmd} {bf_run_prg} {reads_2_file} {vec_name}.bf {k} {hit} {bf_result_2_file} 2> {log_file}")
time_2, mem_2 = parse_log(log_file)

predicted_labels_1 = np.loadtxt(bf_result_1_file, dtype=int)
predicted_labels_2 = np.loadtxt(bf_result_2_file, dtype=int)
tp, tn, fp, fn = calculate_results(true_labels, predicted_labels_1, predicted_labels_2)
precision, recall, f1, accuracy = calculate_metrics(tp, tn, fp, fn)

with open(all_results_file, "w") as result_file:
	result_file.write("READS.1\n")
	result_file.write(f"\ttime:\t{time_1}\n")
	result_file.write(f"\tmem:\t{mem_1}\n")

	result_file.write("READS.2\n")
	result_file.write(f"\ttime:\t{time_2}\n")
	result_file.write(f"\tmem:\t{mem_2}\n")

	result_file.write("CLASS\n")
	result_file.write(f"\ttp:\t{tp}\n")
	result_file.write(f"\ttn:\t{tn}\n")
	result_file.write(f"\tfp:\t{fp}\n")
	result_file.write(f"\tfn:\t{fn}\n")

	result_file.write("METRICS\n")
	result_file.write(f"\tprecision:\t{precision}\n")
	result_file.write(f"\trecall:\t{recall}\n")
	result_file.write(f"\tf1:\t{f1}\n")
	result_file.write(f"\taccuracy:\t{accuracy}\n")
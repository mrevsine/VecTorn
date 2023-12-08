import os
import subprocess
import numpy as np
import re

def calculate_metrics(true_labels, predicted_labels):
    tp = np.sum((true_labels == 1) & (predicted_labels == 1))
    tn = np.sum((true_labels == 0) & (predicted_labels == 0))
    fp = np.sum((true_labels == 0) & (predicted_labels == 1))
    fn = np.sum((true_labels == 1) & (predicted_labels == 0))

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

ref_name = "data/ref/chr22"
vec_name = "data/univec/univec"

bf_build_src = "src/bf/bf_build.cpp"
bf_run_src = "src/bf/bf_run.cpp"
bf_exe_folder = "build"
bf_build_prg = f"{bf_exe_folder}/bf_build"
bf_run_prg = f"{bf_exe_folder}/bf_run"

sim_path = "data/full_sim"
sim_reads_file = f"{sim_path}/reads.1.fastq"
subsample_read_file = f"{sim_path}/subsample.1.fastq"
subsample_class = f"{sim_path}/subsample.class"
subsample_lines = f"{sim_path}/subsample_lines.txt"

bf_result_file = "results/bf/bf_class.txt"
log_file = "results/log"
experiment_results = "results/bf/bf_param_experiments.txt"

if not os.path.isfile(f"{subsample_read_file}"):
	# run_cmd("awk 'NR==FNR {lines[$1]; next} FNR in lines {print; for (i=1; i<=3; i++) {getline; print}}' " + f"{subsample_lines} {sim_reads_file} > {subsample_read_file}")
	print(f"Need to have {subsample_read_file} to run")
	exit
if not os.path.isfile(f"{subsample_class}"):
	run_cmd("awk '/^@/ {print /^@chr22/ ? 0 : 1}' " + f"{subsample_read_file} > {subsample_class}")

if not os.path.exists(bf_exe_folder):
	os.makedirs(f"{bf_exe_folder}/")
if not os.path.isfile(bf_build_prg):
	run_cmd(f"g++ {bf_build_src} -o {bf_build_prg}")
if not os.path.isfile(bf_run_prg):
	run_cmd(f"g++ {bf_run_src} -o {bf_run_prg}")

time_cmd = "/usr/bin/time --format=\"Wall Time: %e\nMax Memory: %M\""

k = 21
hash_size = "100M"
num_threads = 16

run_cmd(f"jellyfish count -m {k} -s {hash_size} -t {num_threads} -o {vec_name}.jf {vec_name}.fa")
run_cmd(f"jellyfish dump {vec_name}.jf | sed '/^>/d' > {vec_name}_kmers.txt")
vec_distinct = get_output_cmd(f"wc -l < {vec_name}_kmers.txt")

# run_cmd(f"jellyfish count -m {k} -s {hash_size} -t {num_threads} -o {ref_name}.jf {ref_name}.fa")
# run_cmd(f"jellyfish dump {ref_name}.jf | sed '/^>/d' > {ref_name}_kmers.txt")
# ref_distinct = get_output_cmd(f"wc -l < {ref_name}_kmers.txt")

# prob_error_ref = 0.001
# run_cmd(f"{bf_build_prg} {ref_name}_kmers.txt {ref_name}.bf {ref_distinct} {prob_error_ref}")

true_labels = np.loadtxt(subsample_class, dtype=int)

with open(experiment_results, "w") as fp:
	fp.write("use_ref\tp\thit\taccuracy\tprecision\trecall\tf1\ttime\tmem\tbf_size\n")

	prob_errors_vec = [1e-2, 1e-4, 1e-6, 1e-8]
	hit_threshold = [2, 4, 6, 8, 10]
	for p in prob_errors_vec:
		run_cmd(f"{bf_build_prg} {vec_name}_kmers.txt {vec_name}.bf {vec_distinct} {p}")
		for hit in hit_threshold:
			# Run without reference bloom filter
			run_cmd(f"{time_cmd} {bf_run_prg} {subsample_read_file} {vec_name}.bf {k} {hit} {bf_result_file} 2> {log_file}")
			time, mem = parse_log(log_file)

			predicted_labels = np.loadtxt(bf_result_file, dtype=int)
			precision, recall, f1, accuracy = calculate_metrics(true_labels, predicted_labels)
			bf_size = get_output_cmd(f"wc -c < {vec_name}.bf")

			fp.write(f"false\t{p}\t{hit}\t{accuracy}\t{precision}\t{recall}\t{f1}\t{time}\t{mem}\t{bf_size}\n")


			# Run with reference bloom filter
			# run_cmd(f"{time_cmd} -p {bf_run_prg} {subsample_read_file} {vec_name}.bf {k} {hit} {bf_result_file} {ref_name}.bf > {log_file}")
			# time, mem = parse_log(log_file)

			# predicted_labels = np.loadtxt(bf_result_file, dtype=int)
			# precision, recall, f1, accuracy = calculate_metrics(true_labels, predicted_labels)
			# bf_size = get_output_cmd(f"wc -c < {vec_name}.bf")

			# fp.write(f"false\t{p}\t{hit}\t{accuracy}\t{precision}\t{recall}\t{f1}\t{time}\t{mem}\t{bf_size}\n")
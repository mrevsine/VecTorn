import os
import subprocess
import re

def run_cmd (cmd):
	print(cmd)
	os.system(cmd)

def get_output_cmd (cmd):
	print(cmd)
	return subprocess.getoutput(cmd)

ref_name = "data/ref/chr22"
vec_name = "data/univec/univec"

results_file = "results/kmer_compare.txt"
# if not os.path.exists("results/"):
# 	os.makedirs("results/")

kmers = [11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
hash_size = "100M"

with open(results_file, 'w') as fp:
	fp.write("k\t#ref\t#vec\t#both\n")
	for k in kmers:
		run_cmd(f"jellyfish count -m {k} -s {hash_size} -o {ref_name}.jf {ref_name}.fa")
		output = get_output_cmd(f"jellyfish stats {ref_name}.jf")
		ref_distinct = int(re.search(r"Distinct:\s+(\d+)", output).group(1))
		run_cmd(f"jellyfish dump {ref_name}.jf | sed '/^>/d' > {ref_name}_kmers.txt")

		run_cmd(f"jellyfish count -m {k} -s {hash_size} -o {vec_name}.jf {vec_name}.fa")
		output = get_output_cmd(f"jellyfish stats {vec_name}.jf")
		vec_distinct = int(re.search(r"Distinct:\s+(\d+)", output).group(1))
		run_cmd(f"jellyfish dump {vec_name}.jf | sed '/^>/d' > {vec_name}_kmers.txt")

		num_both = get_output_cmd(f'bash -c "comm -12 <(sort {ref_name}_kmers.txt) <(sort {vec_name}_kmers.txt) | wc -l"')

		fp.write(f"{k}\t{ref_distinct}\t{vec_distinct}\t{num_both}\n")
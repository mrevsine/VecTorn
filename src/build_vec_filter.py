import os
import subprocess

def run_cmd (cmd):
	print(cmd)
	os.system(cmd)

def get_output_cmd (cmd):
	print(cmd)
	return subprocess.getoutput(cmd)

ref_name = "data/ref/chr22"
vec_name = "data/univec/univec"
bf_build_src = "src/bf_build.cpp"
bf_build_path = "build"
bf_build_prg = f"{bf_build_path}/bf_build"

if not os.path.exists(bf_build_path):
	os.makedirs(f"{bf_build_path}/")

if not os.path.isfile(bf_build_prg):
	os.system(f"g++ {bf_build_src} -o {bf_build_prg}")

k = 21
hash_size = "100M"
num_threads = 16

run_cmd(f"jellyfish count -m {k} -s {hash_size} -t {num_threads} -o {vec_name}.jf {vec_name}.fa")
run_cmd(f"jellyfish dump {vec_name}.jf | sed '/^>/d' > {vec_name}_kmers.txt")
vec_distinct = get_output_cmd(f"wc -l < {vec_name}_kmers.txt")

prob_error_rate = 0.0001
run_cmd(f"{bf_build_prg} {vec_name}_kmers.txt {vec_name}.bf {vec_distinct} {prob_error_rate}")
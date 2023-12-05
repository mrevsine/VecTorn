import os

def __main__():
	# Specify reference seq 
	reference_seq_fa_file = "data/ref/chr22_23M_24M.fa"

	# Specify simulated reads
	simread_file_prefix = "data/full_sim2/reads"
	simread_fastq1 = simread_file_prefix + ".1.fastq"
	simread_fastq2 = simread_file_prefix + ".2.fastq"

	# Bowtie params
	bowtie_file_prefix = "data/ref/chr22_23M_24M"
	bowtie_index_file_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
	bam_file = f"data/ref/"
	
	missing_files = [file for file in bowtie_index_file_suffixes if not os.path.exists(bowtie_file_prefix + file)]
	if missing_files:
		print(f"Building bowtie2 index files for {reference_seq_fa_file}")
		os.system("bowtie2-build " + reference_seq_fa_file + " " + bowtie_file_prefix)
		print("Index built")

	print(f"Aligning {simread_file_prefix} to {reference_seq_fa_file} with bowtie2")
	os.system(f"bowtie2 -x {bowtie_file_prefix} -1 {simread_fastq1} -2 {simread_fastq2} --very-fast -p 16 | samtools view -bS - > out.bam")

if __name__ == "__main__":
	__main__()
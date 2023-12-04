#!/usr/bin/env bash

sim_path=data/full_sim_0error
fq1=${sim_path}/reads.1.fastq
fq2=${sim_path}/reads.2.fastq
ref_name=chr22_23M_24M
ref_fa=data/ref/${ref_name}.fa
ref_index_folder=data/ref/${ref_name}
ref_index=${ref_index}/${ref_name}

if [ ! -d $ref_index_folder ]; then
  bowtie2-build $ref_fa $ref_index -q
fi

n_threads=16

(bowtie2 -x $ref_index -q -1 $fq1 -2 $fq2 --threads $n_threads | samtools view -b - 1>$sim_path/alignment.bam) 2>$sim_path/alignment_log.txt

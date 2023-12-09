#!/usr/bin/env bash

index_dir=data/ref/bwaidx
alignment_dir=data/full_sim
out_dir=results/comparisons/Vecuum

# sort bam file (required by Vecuum)
/usr/local/samtools sort $alignment_dir/alignment.bam -o $alignment_dir/alignment.sorted.bam

# index bam file
/usr/local/samtools index $alignment_dir/alignment.sorted.bam

# run Vecuum
java -jar /usr/local/Vecuum-1.0.1/Vecuum.jar -r $index_dir/chr22_23M_24M.fa -b $alignment_dir/alignment.sorted.bam -B /usr/local/bwa -S /usr/local/samtools -o $out_dir

# extract read names
cut -f 2 $out_dir/*.contaminated.reads | sort | uniq > $out_dir/reads.contaminated.txt
echo "done extracting contaminated reads"

# extract non-contaminated reads
grep "^>" $alignment_dir/reads.1.fa | awk '{sub(/^>/, ""); sub(/\/1$/, ""); print}' | sort > $out_dir/sorted_all_reads.txt
comm -23 $out_dir/sorted_all_reads.txt $out_dir/reads.contaminated.txt > $out_dir/reads.noncontaminated.txt
echo "done extracting non-contaminated reads"

#!/usr/bin/env bash

# list of all paths
vector_dir=data/univec/Univec_Core
vector_db=$vector_dir/db
query_out=results/comparisons/VecScreen/
ref_dir=data/full_sim

# generate vector database if not yet generated
if [[ ! -s $vector_db ]]; then
    /usr/local/ncbi/blast/bin/makeblastdb -in $vector_dir/Univec_Core.fa -dbtype nucl -out $vector_db
fi

# iterate all generated fastq files from the simulation step
for f in $fastq_dir/*.fastq; do 
    # get filename without extension
    filename=$(basename $f)
    filename=${filename%.fastq}
    echo "processing $filename"

    # convert fastq to fasta
    cat $f | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $ref_dir/$filename.fa
    echo "done converting fastq to fasta"

    # query indexed vector database
    /usr/local/ncbi/blast/bin/blastn -query $ref_dir/$filename.fa -db $vector_db -out $query_out/$filename.output.txt -outfmt 6 -num_threads 12
    echo "done querying"

    # extract unique contaminated reads
    cut -f 1 $query_out/$filename.output.txt | sort | uniq > $query_out/$filename.contaminated.txt

    # extract non-contaminated reads
    grep "^>" $ref_dir/$filename.fa | awk 'sub(/^>/, "")' | sort > $query_out/sorted.$filename.txt
    comm -23 $query_out/sorted.$filename.txt $query_out/$filename.contaminated.txt > $query_out/$filename.noncontaminated.txt
    echo "done extracting non-contaminated reads"
done

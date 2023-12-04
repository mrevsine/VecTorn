# VecTorn
Program to identify vector contamination in genomic reads

## Methods

### Creation of a synthetic dataset
#### Recombinant sequeunce simulation
10 genes were identified from human chr22:23M-24M of the GRCH38 release: BCR, DDT, DRICH1, GNAZ, GSTTB2, IGLL1, RAB36, SLC2A11, SMARCB1, and ZNF70. For each, we downloaded their mRNA sequence from NCBI, then simulated converting it to its corresponding CDNA sequence by removing any polyA tails. 10 vectors with an annotated Multiple Cloning Site (MCS) were chosen from the UniVec database: AF113968.2, D63844.1, DQ845760.1, EU186085.1, GU733833.1, JN581376.1, JN900237.1, JX025642.1, U13872.1, and X65329.1. We simulated creating a recombinant sequence from each pair of 10 genes and 10 vectors by inserting the gene's CDNA sequence into the middle of the vector's MCS. For each recombinant sequence we then simulated creating a single SNP within its insert gene region 2 times, for a total of 200 artificial recombinant vectors, each with a single unique mutation.
#### Reference mosaic sequence simulation
10 random loci were selected from human chr22:23M-24M, and a probability distribution of SNP bases was randomly established at each. We then created 20 alternative versions of chr22:23M-24M using simuG by randomly choosing an alternate base at each locus under its corresponding probability table.
#### Paired-end reads simulation
For both the recombinant and reference mosaic sequences, we ran wgsim to randomly simulate paired-end reads of length 150bp with an end-to-end distance of 500 bases. We generated reads to reach a coverage of 250 for our reference mosaic sequences and a coverage of 1000 for our recombinant sequences. The paired-end reads for both sets of sequences were concatenated into one set of paired fastq files.
#### Alignment to the reference region
We aligned the paired fastq files to human chr22:23M-24M using Bowtie2.

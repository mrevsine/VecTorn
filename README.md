# VecTorn
Program to identify vector contamination in genomic reads

## Methods

### Data simulation
10 genes were identified from human Chr22:23M-24M of the GRCH38 release: BCR, DDT, DRICH1, GNAZ, GSTTB2, IGLL1, RAB36, SLC2A11, SMARCB1, and ZNF70. For each, we downloaded their mRNA sequence from NCBI, then simulated converting it to its corresponding CDNA sequence by removing any polyA tails. 10 vectors with an annotated Multiple Cloning Site (MCS) were chosen from the UniVec database: AF113968.2, D63844.1, DQ845760.1, EU186085.1, GU733833.1, JN581376.1, JN900237.1, JX025642.1, U13872.1, and X65329.1. We simulated creating a recombinant sequence from each pair of 10 genes and 10 vectors by inserting the gene's CDNA sequence into the middle of the vector's MCS. For each recombinant sequence we then simulated creating a single SNP within its insert gene region 2 times, for a total of 200 artificial recombinant vectors, each with a single unique mutation.

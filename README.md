# Goals
 

## Introduction

Create a (semi)automated pipeline/method that allows for the detection of genomic features with respect to WGD events and/or transitions from sexual to asexualality.

The following test sets were manually annotated and will serve as positive controls in the automated verison.

#### Test Case 1: 

Dectecting mutation accumulation in a recent WGD event that results in the pseudogenization of duplicated genes.

*Potamopyrgus estuarinus* - pre-duplication state
*Potamopyrgus antipodarum* - closely related species that went through a recent genome duplication event. 

Prior research: 

1. [Angie's Thesis](/files/projects/Codon-Project/articles/honorsthesisakalwies.pdf)


##### Strategy

1. Run BUSCO on transcriptome of *P. estuarinus* and *P. antipodarum*
2. Determine genes found in a single copy and duplicated copies, respectively
3. Blast genome to determine genomic information of genes
4. Blast/BWA raw reads to genomes 
5. de novo assemble raw reads
6. MSA of dupicate pair and single copy
7. Analysis


#### Test Case 2:

Detecting mutation accumulation in a recent (10,000 - 1MYA) transition to asexualality in *Diachasma M.*

*Diachasma muliebre* - asexual




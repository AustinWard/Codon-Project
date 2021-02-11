## Introduction

Create a (semi)automated pipeline that allows for the detection of specific genomic signatures resulting from WGD events and/or transitions from sexual to asexual reproduction. In theory, WGD results in the duplication of the entire gene set, and with that, mutations can accumulate in one copy and avoid selection altogether. This can result in many different outcomes of the duplicated gene (i.e pseudogenization, neo-functionalization and sub-functionalization) that can be detected with comparison to a closely related species in a supposed pre-duplication state. As with sexual vs. asexual, the asexual species has the inability to remove disadvantagous mutations in the absence of sex and recombination. This can result in similar signatures that have been detected in a small set of genes.

The following test sets were manually annotated and will serve as positive controls for the automated verison.

### **Test Case 1:**

Dectecting mutation accumulation in a recent WGD event that results in the pseudogenization of duplicated genes.

#### **Species**

1. *Potamopyrgus estuarinus* - pre-duplication state
2. *Potamopyrgus antipodarum* - closely related species that went through a recent genome duplication event.

#### **Prior research**

1. [Angie's Thesis](/sources/articles/honorsthesis_akalwies.pdf)

![table1](/sources/pictures/sampleout_put001.png)

#### **Strategy**

1. Run BUSCO on transcriptome of *P. estuarinus* and *P. antipodarum*
2. Determine genes found in a single copy and duplicated copies, respectively
3. Blast genome to determine genomic information of genes
4. Blast/BWA raw reads to genomes 
5. de novo assemble raw reads
6. MSA of dupicate pair and single copy
7. Analysis

### **Test Case 2:**

Detecting mutation accumulation in a recent (10,000 - 1MYA) transition to asexualality in *Diachasma M.*

#### **Species**

1. *Diachasma muliebre* - asexual
2. *Diachasma ferrugineum* - sexual
3. *Diachasma alloeum* - outgroup

#### **Prior research**

1. [ET1](/sources/articles/GenomeEvolutionAsexualwasp.pdf)

#### **Strategy**

Similar to that of the first test case. The analysis step is what is being worked on now.

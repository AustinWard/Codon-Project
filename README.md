# Genome Evolution in a Putatively Asexual Wasp (2020)

Eric S. Tvedte et al., Genome evolution in a putatively asexual wasp. BioRxIV https://www.biorxiv.org/content/10.1101/2020.12.23.424202v1.full

## Introduction

Sexual reproduction carries with it many costly attributes (e.g finding a mate, cost of males etc.), however, the long-term survival/adaptations that can be maintained/expunged/generated from recombination alone make up for these costly features. On the other hand, asexual reproducing oraganisms bypass these obstacles at the cost of only being able to pass on traits (mutations) that they themselves accumulate. In theory, evidence for mutation accumlation in asexuals should be detectable at the genomic level with respect to a closely related sexual relative. Evidence for mutations accumlating faster in a asexual species has been shown using *Diachasma mulibre* (asexual) and *D. ferrugineum* (sexual) (Tvedte 2020) and similar findings will be replicated in this project. The goal is to create an automated version of the methods used in this paper so that taxa that meet this asexual/sexual relationship can be analyzed quickly and accurately with the hopes that generalizations can be made with respect to sex and evolution. The results in the table bleow will be replicated to serve as a postive control for the algorithm that generates these counts.

#### **Prior Research Results**

[ET1](/sources/articles/GenomeEvolutionAsexualwasp.pdf)

![table2](/sources/pictures/sample_output002.png)

For the most part, mutations at positions 1 and 2 within a codon result in a non synonymous (dN) mutation - a mutation that changes the amino acid it encodes for. Position 3 mutations (w/ the exception Methionine) result in a synonmous mutation (dS) - a mutation that doesn't result in an amino acid change. The ratio between dN/dS ((1+2)/3)can be used as a marker for evolutionary change, and in this paper, it is hypothesized that it should be greater in *D. mulibre*. 

## **Materials & Methods**

#### Data

The data on the repo (~600 BUSCO genes) was generated differently from the data used to create the table above. The 3,127 genes (nucleotide sequences) were provided by the authors directly and will be used to replicated results similar in the paper. This will serve as way to test the algorithm that produces the counts automatically. From the BUSCO gene set, anotations were collected from running BUSCO --version 4.0.5 on each of the assemblies using metazoa_odb10 db. The shared genes directory contains all COMPLETE single copy busco genes that are shared between all species of wasp. For all species, each gene is it's own file in fasta format. The aligned_proteins directory contains the multiple sequence alignments for each of the fasta files in the proteins directory. These alignements can then be used with the unaligned nucleotide files for the PAL2NAL program. PAL2NAL, which is also used in the paper, allows for a codon aware alignment of the nucleotide sequences and it's output is stored in the codon_align directory. These files are input to a scripts that I'm currently working on to make the counts in the table above.

#### Species

1. *Diachasma muliebre* - asexual - unreleased
2. *Diachasma ferrugineum* - sexual - unreleased
3. *Diachasma alloeum* - outgroup - GCA_001412515.1

#### Methods

* Methods in bold are uinque to using the BUSCO data

**1.** Run [BUSCO](https://gitlab.com/ezlab/busco) on genomes of *D. muliebre*,*D. ferrugineum*, and *D. alloeum* (equalivant to genome annotation)

2. Determine genes found in a single copy.

3. Produce multiple sequence alignments of proteins shared between species - using 3rd party software

4. Produce mutiple sequence alignments of the nucleotide data associated with proteins from (3) using PAL2NAL software (a codon aware alignment tool)

5. Count mutations and group by position of mutation within the codon - ** the mutations will be further classified by Non-Syn vs Synonymous **

Mutations that are analyzed must be unique mutation to the species under investigation (*D. mulibre* vs. *D. ferrugineum*). All this means is that a mutation must have occured after the speciaction event and not in the branch connecting the two nor coming from the branch shared with the outgroup (*D. alloem*)

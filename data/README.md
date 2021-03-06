Data within genomes/wasp was collected from running BUSCO --version 4.0.5 on each of the assemblies using metazoa_odb10 db.
The shared genes directory contains all COMPLETE single copy busco genes that are shared between all species. For all species, each gene is it's own file in fasta format.
The aligned_proteins directory contains the multiple sequence alignments for each of the fasta files in the proteins directory. These alignements can then be used with the unaligned nucleotide files for the PAL2NAL program.
PAL2NAL allows for a codon aware alignment of the nucleotide sequences and it's output is stored in the codon_align directory. These files are input to my python script.

The data I'm reproducing is similar in format.

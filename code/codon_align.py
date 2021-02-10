'''
Input:
	1.  paths to directory containing:
			1. full_table.csv (generated from BUSCO)
			2. "single_copy_busco_sequence/" directory containg .faa and .fna files for completed single copy genes (generated from BUSCO)
Output: 
	1. creates new folder ("shared_genes/") that contain concatenated files for each of the genes sequences (.faa and .fna) 
	   for genes that are shared in all genomes.

The files generated from this script can now be used in pal2nal to make codon aware alignments.

'''

from Bio import SeqIO

from os import sys

if __name__ == "__main__":

	#should rework this function in the future - use dataframe library
	def read_full_table(path, seq_list):
		for line in open(path + "full_table.tsv", 'r'):
			if line[0] != "#" and line != "":
				information = line.split("\t")
				if information[1] == "Complete":
					seq_list.append(information[0])

	def read_fasta(path, gene, seq_list):
		for record in SeqIO.parse(path + "busco_sequences/single_copy_busco_sequences/" + gene, "fasta"):
			seq_list.append(record)
		
	genome1_path = sys.argv[1]
	genome2_path = sys.argv[2]
	genome3_path = sys.argv[3]

	output_path = sys.argv[4]

	genome1_sequences = []
	genome2_sequences = []
	genome3_sequences = []

	read_full_table(genome1_path, genome1_sequences)
	read_full_table(genome2_path, genome2_sequences)
	read_full_table(genome3_path, genome3_sequences)

	print(len(genome1_sequences))
	print(len(genome2_sequences))
	print(len(genome3_sequences))


	gene_intersection = list(set(genome1_sequences) & set(genome2_sequences) & set(genome3_sequences))
	print(len(gene_intersection))


	print(gene_intersection[0:9])

	for gene in gene_intersection:
		nuc_sequeces = []
		aa_sequences = []

		read_fasta(genome1_path, gene + ".faa", aa_sequences)
		read_fasta(genome1_path, gene + ".fna", nuc_sequeces)

		read_fasta(genome2_path, gene + ".faa", aa_sequences)
		read_fasta(genome2_path, gene + ".fna", nuc_sequeces)

		read_fasta(genome3_path, gene + ".faa", aa_sequences)
		read_fasta(genome3_path, gene + ".fna", nuc_sequeces)

		#need to rename header information to include species information
		with open(output_path + "nucleotides/" + gene + "_triplet.fna", "w") as f:
			SeqIO.write(nuc_sequeces, f, "fasta")

		with open(output_path + "proteins/" + gene + "_triplet.faa", "w") as f:
			SeqIO.write(aa_sequences, f, "fasta")

#run output files into pal2nal
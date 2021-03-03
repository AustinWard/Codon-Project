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
from Species import Species

from codon_align import intialize_taxon
from os import sys

import argparse
import os
import json

def read_fasta(path, gene, seq_list):
	for record in SeqIO.parse(path + "/single_copy_busco_sequences/" + gene, "fasta"):
		seq_list.append(record)

def write_fasta(output_path, seq_type, gene, sequences, suffix):
	with open(output_path + seq_type + gene + suffix, "w") as f:
		SeqIO.write(sequences, f, "fasta")

def busco_intersect(config, taxon):

	ortholog_ids = taxon.get_orthologs()

	for gene in ortholog_ids:
		nuc_sequeces = []
		aa_sequences = []
		for s in taxon.species:
			read_fasta(s.busco_path, gene + ".faa", aa_sequences)
			read_fasta(s.busco_path, gene + ".fna", nuc_sequeces)
		
		write_fasta(config['busco_out'], "protein/", gene, aa_sequences, "_triplet.faa")
		write_fasta(config['busco_out'], "nucleotide/", gene, nuc_sequeces, "_triplet.fna")

def main(args):

	if os.path.exists(args.config):
		with open(args.config) as f:
			config = json.load(f)
	else: 
		print("Config file not found")
		sys.exit()

	taxon = intialize_taxon(config)
	busco_intersect(config, taxon)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Counting positional mutations within codons')
	parser.add_argument("config", help='Configuration file for run. See example config for details.')

	main(parser.parse_args())
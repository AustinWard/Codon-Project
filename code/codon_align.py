from Bio import AlignIO, SeqIO
from CodonTable import CodonTable
from Species import Species
from Taxon import Taxon

import argparse
import sys
import json
import os

def intialize_taxon(config):

	#create species
	species = []
	for s in config['species']:
		species.append(Species(s['name'], config['taxon'], s['busco_data'], s['outgroup']))

	#create taxon
	taxon = Taxon(config['taxon'], species)

	return(taxon) 

def count_mutations(taxon):

		outlier = None
		in_group = []

		for s in taxon.species:
			if s.outlier:
				outlier = s
			else: in_group.append(s)
		
		if outlier is not None:
			#cycling through codons and comparing them to the codons of other species within the taxon
			for out_codon in outlier.codons:

				taxon.total_length+=3
				current_codons = []

				found_match = False #if one of the in-group amino acides matches the outlier - we need at least one
	
				if out_codon.aa != '-':
					found_match = True
					for seq in in_group:
						#TO DO: should include test as to whether the lengths of the protein and nucleotide sequences are of the same length
						current_aa = seq.codons[out_codon.pos].aa
						#print(seq.codons[out_codon.pos])
						
						if current_aa == '-':
							found_match = False
							break

						current_codons.append(seq.codons[out_codon.pos])

				if found_match:

					analyze_positions = [False, False, False]
					c1n1 = current_codons[0].n1
					c1n2 = current_codons[0].n2
					c1n3 = current_codons[0].n3

					counter = 1
					while counter < len(current_codons):
						if c1n1 != current_codons[counter].n1:
							analyze_positions[0] = True
						if c1n2 != current_codons[counter].n2:
							analyze_positions[1] = True
						if c1n3 != current_codons[counter].n3:
							analyze_positions[2] = True
						counter+=1
		
					taxon.counted_positions+=3

					for codon in current_codons:
						#Looking for Syn. mutation - all mutations have to be syn. by definition
						found_syn = False
						found_nsyn = False

						if out_codon.aa == codon.aa:
							if codon.n1 != out_codon.n1 and analyze_positions[0]:
								found_syn = True
								codon.species.sp1+=1
							if codon.n2 != out_codon.n2 and analyze_positions[1]:
								found_syn = True
								codon.species.sp2+=1
							if codon.n3 != out_codon.n3 and analyze_positions[2]:
								found_syn = True
								codon.species.sp3+=1
							if found_syn:
								codon.species.syn_count+=1

						#Non-Syn mutations in the 3rd position - example 
						elif((codon.n1 == out_codon.n1) and (codon.n2 == out_codon.n2)) and analyze_positions[2]:
							found_nsyn = True
							codon.species.nsp3+=1

						#Look for Non. Syn mutations - not all mutations are non-syn. - checks for these 
						else:

							#look at first two positions - any mutation in the first two positions will
							#always result in a non-syn mutation
							if codon.n1 != out_codon.n1 and analyze_positions[0]:
								found_nsyn = True
								codon.species.nsp1+=1
							if codon.n2 != out_codon.n2 and analyze_positions[1]:
								codon.species.nsp2+=1
								found_nsyn = True
							if codon.n3 != out_codon.n3 and analyze_positions[2]:
								found_nsyn = True
								if out_codon.aa == 'M':
									if analyze_positions[2]:
										codon.species.nsp3+=1
								if out_codon.aa == 'I':
									if codon.n3 == 'G':
										if analyze_positions[2]:
											codon.species.nsp3+=1
									else: 
										if analyze_positions[2]:
											codon.species.sp3+=1

								#aa that are encoded by codons that can have anything at the 3rd position - replace next line with regex
								elif ((out_codon.n2 == 'C') or ((out_codon.n1 in ['C', 'G']) and (out_codon.n2 in ['T', 'G']))):
									if analyze_positions[2]:
										codon.species.sp3+=1

								#aa that can have T/C at 3rd position OR A/G at 3rd position
								else:
									if (out_codon.n3 in ['T', 'C']):
										if (codon.n3 in ['T', 'C']):
											if analyze_positions[2]:
												codon.species.sp3+=1
										else:
											if analyze_positions[2]:
												codon.species.nsp3+=1
									#A or G in 3rd position
									else:
										if (codon.n3 in ['A', 'G']):
											if analyze_positions[2]:
												codon.species.sp3+=1
										else:
											if analyze_positions[2]:
												codon.species.nsp3+=1

						if found_nsyn: codon.species.nsyn_count+=1

def read_seqs(config, taxon, codon_table):

	#looping through genes to be analyzed if directory exists
	if os.path.isdir(config['input_directory']):
		count = 0
		genes = os.listdir(config['input_directory'])
		taxon.total_genes = len(genes)

		while count < taxon.total_genes:
			file_path = os.path.join(config['input_directory'], genes[count])
			with open(file_path) as f:
				nuc_msa = AlignIO.read(f, 'fasta')

			if len(nuc_msa) < config['consider']:
				print("The number of sequences in " + genes[count] + " is less than what is expected.")
				sys.exit()

			seq_count = 0
			while seq_count < config['consider']:
				taxon.species[seq_count].nucleotide = nuc_msa[seq_count]
				taxon.species[seq_count].add_codons(codon_table)

				seq_count+=1

			count_mutations(taxon)

			count+=1

def busco_intersect(taxon):
	pass


def main(args):

	codon_table = CodonTable() #will allow user to select a specific codon table in the future

	## read json if it exists
	if os.path.exists(args.config):
		with open(args.config) as f:
			config = json.load(f)
	else: 
		print("Config file not found")
		sys.exit()

	taxon = intialize_taxon(config)
	read_seqs(config, taxon, codon_table)

	for s in taxon.species:
		print(s)
	print(taxon)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Counting positional mutations within codons')

	parser.add_argument("config", help='Configuration file for run. See example config for details.')
	#parser.add_argument('--p_aligns', nargs='+', help='Protein alignment file in fasta format (can add more formats later).')
	#parser.add_argument('--c_aligns', nargs='+', required=True, help='Corresponding codon alignment file from pal2nal.')
	
	main(parser.parse_args())

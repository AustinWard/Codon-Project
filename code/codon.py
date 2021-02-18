import argparse
import sys
import json
import os

from Bio import AlignIO
from CodonTable import CodonTable 

class Codon():

	def __init__(self, n1, n2, n3, pos, species, aa='None'):
		self.aa = aa
		self.n1 = n1
		self.n2 = n2
		self.n3 = n3
		self.pos = pos

		self.species = species

	def __str__(self):
		return("Species: " + str(self.species.name) + " Position: " + str(self.pos) + " Amino Acid: " + self.aa + "\nNucleotides: " + self.n1 + self.n2 + self.n3)

	def add_aa(self, codon_table):
		codon = self.n1 + self.n2 + self.n3
		codon = codon.upper()
		self.aa = codon_table.get_aa(codon)

class Species():

	def __init__(self, name, taxa, outlier=False,):
		self.name = name
		self.taxa = taxa
		self.outlier = outlier
		
		self.protein = None
		self.nucleotide = None
		self.codons = []

		#total synonymous mutations 1,2,3rd position of a codon
		self.sp1 = 0
		self.sp2 = 0
		self.sp3 = 0
		self.syn_count = 0

		#total Nonsynonymous mutations 1,2,3rd position of a codon
		self.nsp1 = 0
		self.nsp2 = 0
		self.nsp3 = 0
		self.nsyn_count = 0

	def add_codons(self, codon_table):

		self.codons = []

		seq_length = len(self.nucleotide)

		if seq_length % 3 != 0:
			print("sequence isn't a multiple of 3")
			sys.exit()
		else:
			codon_count = seq_length / 3
			index = 0
			nuc_start = 0
			while index < codon_count:

				n1 = self.nucleotide[nuc_start]
				n2 = self.nucleotide[nuc_start+1]
				n3 = self.nucleotide[nuc_start+2]

				codon = Codon(n1, n2, n3, index, self)
				try:
					codon.add_aa(codon_table)
				except:
					codon.aa = '-'
				
				self.codons.append(codon)

				nuc_start+=3
				index+=1

	def __str__(self):
		return(self.taxa + " " +self.name + "\n" +
			"_____________________________\n" +
			"_____________________________\n" +
			"Nonsynonymous: \n" +
			"P1: " + str(self.nsp1) + "\n" +
			"P2: " + str(self.nsp2) + "\n" +
			"P3: " + str(self.nsp3) + "\n\n" +
			"Nonsynonymous Count: " + str(self.nsyn_count) + "\n"
			"_____________________________\n" +
			"Synonymous: \n" +
			"P1: " + str(self.sp1) + "\n" +
			"P2: " + str(self.sp2) + "\n" +
			"P3: " + str(self.sp3) + "\n\n" +
			"Synonymous Count: " + str(self.syn_count) + "\n" +
			"_____________________________\n" +
			"_____________________________\n")

class Group():

	def __init__(self, name, species):
		self.name = name 						#name of taxa 
		self.species = species 					#list of species objects

		self.total_genes = 0
		self.total_length = 0
		self.counted_positions = 0				#number of nucleotide positions looked at

	def __str__(self):
		return("Group: " + self.name + "\n" + 
			"Genes: " + str(self.total_genes) + "\n" +
			"Total Length: " + str(self.total_length) + "\n"
			"Counted Positions: " + str(self.counted_positions) + "\n")

	def count_mutations(self):

		outlier = None
		in_group = []

		for s in self.species:
			if s.outlier:
				outlier = s
			else: in_group.append(s)
		
		if outlier is not None:
			#cycling through codons and comparing them to the codons of other species within the taxa
			for out_codon in outlier.codons:

				self.total_length+=3
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
		
					self.counted_positions+=3

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

def main(args):

	codon_table = CodonTable()

	#parse configuration file

	## read json if it exists
	if os.path.exists(args.config):
		with open(args.config) as f:
			config = json.load(f)

	else: 
		print("Config file not found")
		sys.exit()

	#create species
	species = []
	for s in config['species']:
		species.append(Species(s['name'], config['group'], s['outgroup']))

	if len(species) > config['consider']:
		print("The number of individuals must be <= to the number of sequences to consider (i.e number of sequences in fasta file(s))")
		sys.exit()

	#create group
	group = Group(config['group'], species)

	#looping through genes to be analyzed if directory exists
	if os.path.isdir(config['input_directory']):
		count = 0
		genes = os.listdir(config['input_directory'])
		group.total_genes = len(genes)

		while count < group.total_genes:
			file_path = os.path.join(config['input_directory'], genes[count])
			with open(file_path) as f:
				nuc_msa = AlignIO.read(f, 'fasta')

			if len(nuc_msa) < config['consider']:
				print("The number of sequences in " + genes[count] + " is less than what is expected.")
				sys.exit()

			seq_count = 0
			while seq_count < config['consider']:
				species[seq_count].nucleotide = nuc_msa[seq_count]
				species[seq_count].add_codons(codon_table)

				seq_count+=1

			group.count_mutations()

			count+=1 
			
	for s in species:
		print(s)
	print(group)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Counting positional mutations within codons')

	parser.add_argument("config", help='Configuration file for run. See example config for details.')
	#parser.add_argument('--p_aligns', nargs='+', help='Protein alignment file in fasta format (can add more formats later).')
	#parser.add_argument('--c_aligns', nargs='+', required=True, help='Corresponding codon alignment file from pal2nal.')
	
	main(parser.parse_args())

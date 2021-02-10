import argparse

from Bio import AlignIO

class Codon():

	def __init__(self, aa, n1, n2, n3, pos, species):
		self.aa = aa
		self.n1 = n1
		self.n2 = n2
		self.n3 = n3
		self.pos = pos

		self.species = species

	def __str__(self):
		return("Species: " + str(self.species.name) + " Position: " + str(self.pos) + " Amino Acid: " + self.aa + "\nNucleotides: " + self.n1 + self.n2 + self.n3)

class Species():

	def __init__(self, name, outlier=False):
		self.name = name
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

	def add_codons(self):

		self.codons = []

		seq_length = len(self.nucleotide)

		#the protein sequences should be exactly a third of length of the codon aligned nucleotide sequence
		if (3*len(self.protein)) == seq_length:
			
			for index, aa in enumerate(self.protein):
				
				nuc_start = index*3

				n1 = self.nucleotide[nuc_start]
				n2 = self.nucleotide[nuc_start+1]
				n3 = self.nucleotide[nuc_start+2]

				codon = Codon(aa, n1, n2, n3, index, self)

				self.codons.append(codon)

	def __str__(self):
		return(self.name + "\n" +
			"_____________________________\n" +
			"Nonsynonymous: \n" +
			"P1: " + str(self.nsp1) + "\n" +
			"P2: " + str(self.nsp2) + "\n" +
			"P3: " + str(self.nsp3) + "\n" +
			"_____________________________\n" +
			"Nonsynonymous Count: " + str(self.nsyn_count) + "\n"
			"Synonymous: \n" +
			"P1: " + str(self.sp1) + "\n" +
			"P2: " + str(self.sp2) + "\n" +
			"P3: " + str(self.sp3) + "\n" +
			"Synonymous Count: " + str(self.syn_count) + "\n" +
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

					self.counted_positions+=3

					for codon in current_codons:
						#Looking for Syn. mutation - all mutations have to be syn. by definition
						found_syn = False
						found_nsyn = False

						if out_codon.aa == codon.aa:
							if codon.n1 != out_codon.n1:
								found_syn = True
								codon.species.sp1+=1
							if codon.n2 != out_codon.n2:
								found_syn = True
								codon.species.sp2+=1
							if codon.n3 != out_codon.n3:
								found_syn = True
								codon.species.sp3+=1
							if found_syn:
								codon.species.syn_count+=1

						#Non-Syn mutations in the 3rd position - example 
						elif((codon.n1 == out_codon.n1) and (codon.n2 == out_codon.n2)):
							found_nsyn = True
							codon.species.nsp3+=1

						#Look for Non. Syn mutations - not all mutations are non-syn. - checks for these 
						else:

							#look at first two positions - any mutation in the first two positions will
							#always result in a non-syn mutation
							if codon.n1 != out_codon.n1:
								found_nsyn = True
								codon.species.nsp1+=1
							if codon.n2 != out_codon.n2:
								codon.species.nsp2+=1
								found_nsyn = True

							if codon.n3 != out_codon.n3:
								found_nsyn = True
								if out_codon.aa == 'M':
									codon.species.nsp3+=1
								if out_codon.aa == 'I':
									if codon.n3 == 'G':
										codon.species.nsp3+=1
									else: codon.species.sp3+=1

								#aa that are encoded by codons that can have anything at the 3rd position - replace next line with regex
								elif ((out_codon.n2 == 'C') or ((out_codon.n1 in ['C', 'G']) and (out_codon.n2 in ['T', 'G']))):
									codon.species.sp3+=1

								#aa that can have T/C at 3rd position OR A/G at 3rd position
								else:
									if (out_codon.n3 in ['T', 'C']):
										if (codon.n3 in ['T', 'C']):
											codon.species.sp3+=1
										else:
											codon.species.nsp3+=1
									#A or G in 3rd position
									else:
										if (codon.n3 in ['A', 'G']):
											codon.species.sp3+=1
										else:
											codon.species.nsp3+=1

						if found_nsyn: codon.species.nsyn_count+=1

def main(args):

	if len(args.p_aligns) == len(args.c_aligns):

		#need to be able to automate this
		species = []
		species.append(Species('DA_WASP', True))
		species.append(Species('DF_WASP'))
		species.append(Species('DM_WASP'))

		#and this
		group = Group('Wasp', species)

		count = 0
		#looping through genes
		while count < len(args.p_aligns):

			p_file = open(args.p_aligns[count])
			c_file = open(args.c_aligns[count])

			p_align = AlignIO.read(p_file, 'fasta')
			codon_align = AlignIO.read(c_file, 'fasta')

			c_file.close()
			p_file.close()

			#same number of sequences and number of sequences is the same as the number of species
			if len(p_align) == len(codon_align) and len(species) == len(p_align):
				
				seq_count = 0
				group.total_genes+=1

				while seq_count < len(p_align):

					species[seq_count].protein = p_align[seq_count]
					species[seq_count].nucleotide = codon_align[seq_count]

					species[seq_count].add_codons()

					seq_count+=1
				
				group.count_mutations()

			count+=1 

	for s in species:
		print(s)
	print(group)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Counting positional mutations within codons')

	parser.add_argument('--p_aligns', nargs='+', required=True, help='Protein alignment file in fasta format (can add more formats later).')
	parser.add_argument('--c_aligns', nargs='+', required=True, help='Corresponding codon alignment file from pal2nal.')
	
	main(parser.parse_args())

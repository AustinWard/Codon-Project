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
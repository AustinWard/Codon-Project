import unittest

from Bio import AlignIO
from codon import Species, Group, Codon
from CodonTable import CodonTable

class TestCounts(unittest.TestCase):

	def setUp(self):

		codon_table = CodonTable()

		self.species = []

		self.species.append(Species('TAXON_1', 'taxa', True))
		self.species.append(Species('TAXON_2', 'taxa'))
		self.species.append(Species('TAXON_3', 'taxa'))

		self.group = Group('Group', self.species)


		p_file = open("../data/test/test_data_triplet.faa")
		c_file = open("../data/test/test_data_triplet.fna")

		p_align = AlignIO.read(p_file, 'fasta')
		codon_align = AlignIO.read(c_file, 'fasta')

		c_file.close()
		p_file.close()

		seq_count = 0

		while seq_count < len(p_align):

			self.species[seq_count].protein = p_align[seq_count]
			self.species[seq_count].nucleotide = codon_align[seq_count]

			self.species[seq_count].add_codons(codon_table)

			seq_count+=1

		self.group.count_mutations()
	#taxon 2 nonsyn p1-3
	def test_s2_nsp1(self):
		self.assertEqual(self.species[1].nsp1, 2)
	def test_s2_nsp2(self):
		self.assertEqual(self.species[1].nsp2, 1)
	def test_s2_nsp3(self):
		self.assertEqual(self.species[1].nsp3, 1)

	#taxon 2 syn. p1-3
	def test_s2_sp1(self):
		self.assertEqual(self.species[1].sp1, 2)
	def test_s2_sp2(self):
		self.assertEqual(self.species[1].sp2, 1)
	def test_s2_sp3(self):
		self.assertEqual(self.species[1].sp3, 3)

	#taxon 3 nonsyn p1-3
	def test_s3_nsp1(self):
		self.assertEqual(self.species[2].nsp1, 2)
	def test_s3_nsp2(self):
		self.assertEqual(self.species[2].nsp2, 3)
	def test_s3_nsp3(self):
		self.assertEqual(self.species[2].nsp3, 0)

	#taxon 3 syn p1-3
	def test_s3_sp1(self):
		self.assertEqual(self.species[2].sp1, 0)
	def test_s3_sp2(self):
		self.assertEqual(self.species[2].sp2, 0)
	def test_s3_sp3(self):
		self.assertEqual(self.species[2].sp3, 4)

if __name__ == "__main__":
	unittest.main()
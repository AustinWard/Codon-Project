from Codon import Codon

class Species():

	def __init__(self, name, taxa, busco_path, outlier=False):
		self.name = name
		self.taxa = taxa
		self.outlier = outlier
		
		self.protein = None
		self.nucleotide = None
		self.codons = []

		self.busco_path = busco_path

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

	def busco_complete(self, path=""):

		if path == "":
			path = self.busco_path
		seq_list = set()
		for line in open(path + "/full_table.tsv", 'r'):
			if line[0] != "#" and line != "":
				information = line.split("\t")
				if information[1] == "Complete":
					seq_list.add(information[0])
		return(seq_list)

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
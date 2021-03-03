from Species import Species

class Taxon():

	def __init__(self, name, species):
		self.name = name 						#name of taxa 
		self.species = species 					#list of species objects

		self.total_genes = 0
		self.total_length = 0
		self.counted_positions = 0				#number of nucleotide positions looked at

	def get_orthologs(self):
		shared_genes = set()
		for s in self.species:
			complete_genes = s.busco_complete()

			#if set is empty - intialize the shared_genes as the first set of genes
			if len(shared_genes) == 0:
				shared_genes = complete_genes
			else:
				shared_genes = shared_genes.intersection(complete_genes)
		return(shared_genes)

	def __str__(self):
		return("Taxon: " + self.name + "\n" + 
			"Genes: " + str(self.total_genes) + "\n" +
			"Total Length: " + str(self.total_length) + "\n"
			"Counted Positions: " + str(self.counted_positions) + "\n")
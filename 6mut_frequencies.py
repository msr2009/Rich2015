"""
6mut_frequencies.py

genotypes paired-end reads for 6mut combinations

Matt Rich, 1/2015
"""

def main(fastq, fastq2, fastq3, positions, exclude):
		
	dat = []
	for read in read_fastq_multi([option.fastq, option.fastq2, option.fastq3]):
		
		#trim and RC reverse read
		revread = reverse_fastq(trim_fastq_length(read[1], 1, 106))
				
		#merge read1 and read2
		merged = zip(*[ merge_reads(read[0][1][x], read[0][2][x], revread[1][x], revread[2][x]) for x in range(len(revread[1])) ])

		#genotype each position
		pos = sorted(positions.keys()) #make a list of the sorted keys
		g = [ merged[0][x] for x in pos ]
		q = [ merged[1][x] for x in pos ]
		
		#check which genotype -- ref, alt, or other	
		genotypes = [ genotype_sites(g[x], positions[pos[x]]) for x in range(len(pos)) ]
		
		#add barcode genotype
		genotypes.append(genotype_sites(read[2][1][5], ["T","G"]))	
		
		if exclude:
			if "N" not in genotypes:
				print ''.join([ str(x) for x in genotypes ])
		else:
			print ''.join([ str(x) for x in genotypes ])

def genotype_sites(base, poss):
	try:
		return poss.index(base)
	except ValueError:
		return "N"

def merge_reads(base1, qual1, base2, qual2):
	if base1 == base2:
		return base1, qual1+qual2
	else:
		if qual1 > qual2:
			return base1, qual1
		else:
			return base2, qual2

if __name__ == "__main__":
	
	from optparse import OptionParser
	from fastq_tools import read_fastq_multi, trim_fastq_length, reverse_fastq

	parser = OptionParser()
	parser.add_option('--fastq', action = 'store', type = 'string', dest = 'fastq', help = "fastq file containing forward reads")
	parser.add_option('--fastq2', action = 'store', type = 'string', dest = 'fastq2', help = "fastq file containing reverse reads")
	parser.add_option('--fastq3', action = 'store', type = 'string', dest = 'fastq3', help = "fastq file containing index reads")
	parser.add_option('--positions', action = 'store', type = 'string', dest = 'positions', help = 'positions to genotype')
	parser.add_option('--print_all', action = 'store_false', dest = 'exclude', help = 'print all genotypes (including Ns?)', default = True)
	(option, args) = parser.parse_args()
	
	positions = {}

	if option.positions == None:	
		positions = {0:["T","A"], 45:["A","G"], 54:["C","T"], 86:["T","C"], 105:["T","G"]}
				
	main(option.fastq, option.fastq2, option.fastq3, positions, option.exclude)

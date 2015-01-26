"""
fimo_newsites.py

using output from FIMO search against motifs, finds positions that are putatively 
creating new sites based on significance of match

Matt Rich, 9/2014
"""

def main(fimo, p, missing):
	# read fimo text file into data structure
	# f[TF;start;stop;strand] = { pos: {"WT":N, "A":pA, "C":pC, "G":pG, "T":pT } }
	f = {}
	for line in open(fimo, "r"):
		l = line.strip().split()	
		if line.startswith("#"):
			print '\t'.join(["#motif", "sequence", "start", "stop", "strand", "p", "p_wt", "p_diff"])
		else:
			pos = l[1][2:-3] 
			ref = l[1][-3]
			obs = l[1][-1]
			motif = ":".join([ l[0], l[2], l[3], l[4], pos ])
#			print pos, ref, obs, motif
			#if we haven't seen the motif match before:
			if motif not in f:
				f[motif] = {}
				#initialize entry in dictionary
				f[motif][pos] = { "WT": ref, "A":missing, "C":missing, "G":missing, "T":missing } 
				#update with the current p-value
				f[motif][pos][obs] = float(l[6])		
			else:
				f[motif][pos][obs] = float(l[6])		

	#calculate p_WT/p_var for all sites
	for m in f:
		m_data = m.split(":")
		for pos in f[m]:
			wt = f[m][pos]["WT"]
			for variant in f[m][pos]:
				if variant != "WT" and variant:
					p_diff = f[m][pos][wt] / f[m][pos][variant]		
					print "\t".join( [ m_data[0], "n." + str(int(pos)+1) + wt + ">" + variant, m_data[1], m_data[2], m_data[3], \
						str(f[m][pos][variant]), str(f[m][pos][wt]), str(p_diff) ] )				

if __name__ == "__main__":
	
	from optparse import OptionParser
	
	parser = OptionParser()

	parser.add_option("--fimo", action = "store", type = "string", dest = "fimo", help = "FIMO output (fimo.txt)")
	parser.add_option("-p", action = "store", type = "float", dest = "pvalue", help = "threshold for significant match", default = 0.0005)
	parser.add_option("--missing", action = "store", type = "float", dest = "missing", help = "p-value for missing positions", default = 0.01)	
	(option, args) = parser.parse_args()

	main(option.fimo, option.pvalue, option.missing)

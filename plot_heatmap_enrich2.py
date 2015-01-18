"""
Plots a heatmap based on the pandas output of Enrich v2.0

The basis for the code to plot the heatmap taken from Dave Young. 
All else (parameterization, barplots, colormap changes, etc..) by MR.

Matt Rich, 7/2014
"""

def main(df, wt_seq, dna, figout, offset, subset, set_min, set_max, xlim):
	
	#set the WT sequence of your region of interest
	WT_sequence=wt_seq
	if dna == False:
		aa = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'R', 'H', 'K', 'D', 'E', '*']
	else:
		aa = ['A', 'C', 'G', 'T']
	
	#This next line gives you the position of each WT amino acid in your "aa" list. 
	#We use this later to mark the WT positions in the heatmap. We reverse all our 
	#lists so that the heatmap is drawn top to bottom (default is bottom to top)
	WT_aa_index = [aa.index(i) for i in WT_sequence] 
	
	#Columns will denote position (first 100 positions here)
	row_labels = aa
	ticks = np.arange(0, len(WT_sequence), 10)	
	
	#make a pivot table of all the data 
	### KINDA HACKY -- THIS ASSUMES THERE'S AT LEAST ONE DATAPOINT FOR EACH POSITION
	data = df.pivot("pos", "mut", "dat")
#	print df.count()

	#also, here we can set the min and max values for our colorbars
	if set_max == None:	
		vmax=max(df["dat"])
		print "vmax: " + str(vmax)
	else:
		vmax=set_max
	
	if set_min == None:
		vmin=min(df["dat"])
		print "vmin: " + str(vmin)
	else:
		vmin=set_min

	column_labels = data.index.values - 1	

	#also make dataframe for mean and sd for each position	
	data2 = pd.DataFrame()
	data2["mu"] = data.mean(1)
	data2["sd"] = data.std(1)
	
	#set heatmap color scheme
	BlueWhiteYellow = mpl.colors.LinearSegmentedColormap.from_list( "BlueWhiteYellow", ((0, "blue"),(.5, "lightgrey"),(1, "yellow")) )
	my_cmap = remappedColorMap(BlueWhiteYellow, start=0, midpoint=abs(vmin)/(vmax+abs(vmin)), stop=1 ) 

	#set color for masked array entries
	my_cmap.set_bad('white') 
	wt_col = "lightblue"
	#make a normalize instance for all the plots
	my_norm = mpl.colors.Normalize(vmin, vmax)
	#make a new Series based on the colormap to define a color for each mean value
	hbar_cols = cm.ScalarMappable(mpl.colors.Normalize(vmin, vmax), my_cmap).to_rgba(data2["mu"])
	
	fig = plt.figure()
	tlen=15
	if subset != None:
		fig.set_size_inches(4,6)
		axarr = plt.subplot()
		data = data.loc[subset[0]:subset[1]+1]
		column_labels = data.index.values-1
#		print data
		tlen=5
		
	else:
		fig.set_size_inches(22,42)
		gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 1]) 
		axarr = plt.subplot(gs[0])
#		print data	
	
	#make heatmap 
	heatmap = axarr.imshow(data, norm=my_norm, interpolation = 'none', aspect='auto', cmap=my_cmap, extent=[0, len(aa), max(column_labels)+1, min(column_labels)])

	#set y-axis
	ytickloc = mpl.ticker.MultipleLocator(10)
	ytickform = mpl.ticker.ScalarFormatter()
	axarr.yaxis.set_major_locator(ytickloc)
	axarr.yaxis.set_major_formatter(ytickform)
	axarr.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
	axarr.set_xticklabels(row_labels, minor=False)
	axarr.set_xticks(np.arange(data.shape[1])+.5, minor=False)
	axarr.tick_params(axis='both', which='major', labelsize=26, length=tlen)
	axarr.tick_params(axis='both', which='both', direction='out')
	axarr.invert_yaxis()
#	axarr.set_frame_on(False)
	axarr.grid(False)
	#Mark the WT positions
	for i in column_labels:
		axarr.add_patch(mpl.patches.Rectangle((WT_aa_index[i],i), 1, 1, fill=True, facecolor=wt_col, edgecolor="none")) 
	
	if subset == None:
		#make mean/sd barplot
		axarr2 = plt.subplot(gs[1], sharey=axarr)
		bars = axarr2.barh(data2.index, data2["mu"], linewidth=0.5, color=hbar_cols, align='center')
		if xlim != None:
			plt.xlim(xlim)
		axarr2.yaxis.set_major_locator(ytickloc)
		axarr2.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
	#	axarr2.errorbar(data2.index, df_mean, yerr = df_sd, capsize = 0, color="r")
		axarr2.set_xlabel("Average wt-normalized slope")
		axarr2.set_ylim(min(column_labels), max(column_labels))
		axarr2.tick_params(axis='both', which='major', labelsize=26, length=15)
		axarr2.invert_xaxis()
#		plt.setp(axarr2.get_yticklabels(), visible=False)
		
	fig.colorbar(heatmap, ticks=[-.2,-.15,-0.1,-0.05,0,0.05,0.1,.15,.2,.25,.3,.35,.4])
#	fig.colorbar(heatmap, ticks=[-0.1,-.075,-0.05,-.025,0,.025,0.05,.075,0.1])
#	fig.title("Score", fontsize=60)
	plt.tight_layout()
	
	plt.savefig(figout)

def preprocess(df, datcol, mincounts, takelog):
	#split HGVS into pos, mut
	new_df = df[df["count.0"] >= mincounts]["sequence"].transpose().apply(func=HGVS_to_PosMut)

	new_df.columns = ["pos", "ref", "mut"]
#	print new_df
	if takelog:
		new_df["dat"] = df[datcol].apply(np.log2)
	else:
		new_df["dat"] = df[datcol]
	return new_df
	
def HGVS_to_PosMut(hgvs):
	if hgvs != "_wt":
		posmut = re.match('n\.(-?\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs).groups()
		return pd.Series([int(posmut[0]), posmut[1], posmut[2]])
	else:
		pass

if __name__ == "__main__":
	
	from optparse import OptionParser
	
	import matplotlib as mpl
	from matplotlib import gridspec
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm

	from matplotlib_tools import remappedColorMap
	
	import pandas as pd
	import numpy as np	
	import re
	
	parser = OptionParser()
	parser.add_option('--data', action = 'store', type = 'string', dest = 'data', help = "tab-delimited dataframe")
	parser.add_option('--dna', action = 'store_true', dest = 'dna', help = "is DNA sequence?", default = False)
	parser.add_option('--log', action = 'store_true', dest = 'log', help = "take log2 of data?", default = False)
	parser.add_option('--min-counts', action = 'store', type = 'int', dest = 'mincounts', help = 'minimum number of input counts', default = 0)
	parser.add_option('--wt', action = 'store', type = 'string', dest = 'wt_seq', help = "wildtype sequence")
	parser.add_option('--out', '-o', action = 'store', type = 'string', dest = 'figout', help = 'output filename', default='fig.png')
	parser.add_option('--datcol', action = 'store', type = 'string', dest = 'datcol', help = 'column to plot')
	parser.add_option('--offset', action = 'store', type = 'int', dest = 'offset', help = 'axis offset', default = 0)
	parser.add_option('--subset', action = 'store', type = 'string', dest = 'subset', help = 'comma-delimited subset of data to plot as heatmap', default = None)
	parser.add_option('--min', action = 'store', type = 'float', dest = 'set_min', help = 'minimum value for colorbar', default = None)
	parser.add_option('--max', action = 'store', type = 'float', dest = 'set_max', help = 'maximum value for colorbar', default = None)
	parser.add_option('--xlim', action = 'store', type = "string", dest = 'xlim', help = 'limits of barplot axis', default=None)
	
	(option, args) = parser.parse_args()
	
	subset = option.subset
	if subset != None:
		subset = [int(x) for x in subset.split(",")]
	
	xlim = option.xlim
	if option.xlim != None:
		xlim = [float(x) for x in option.xlim.split(",")]
	
	main(preprocess(pd.read_csv(option.data, sep="\t", header=0), option.datcol, option.mincounts, option.log), 
		option.wt_seq.upper(), option.dna, option.figout, option.offset, subset, option.set_min, option.set_max, xlim)	



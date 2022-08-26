import os
import math
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import itertools
from matplotlib.lines import Line2D
import matplotlib.font_manager
from matplotlib import rcParams, rcParamsDefault
import numpy as np;
from matplotlib import colors
from cycler import cycler
from math import log

# Global Parameter to enable or disable LaTeX
global_latex_en = True

# Global map: label_in_csv -> label_in_legend
global_label_map = {
					"BASE" 		: "Exact",
					"BF" 		: "ProbGraph (BF)",
					"1H" 		: "ProbGraph (MH)",
					"REDEX" 	: "Reduced Execution",
					"PGP" 		: "Partial Graph Proc.",
		 			"COLORFUL"	: "Colorful",
		 			"DOULION"	: "Doulion",
					"BSP1"		: "AutoApprox1",
					"BSP2"		: "AutoApprox2",
}

# Global map: label_in_csv -> marker
global_marker_map = {
					"BF"		: '^', 
		 			"1H"		: 's', 
		 			"BASE"		: '*', 
					"REDEX"		: 'd',
					"PGP"		: 'D',
		 			"COLORFUL"	: 'o',
		 			"DOULION"	: 'X',
					"BSP1"		: "p",
					"BSP2"		: "h",
}


# Global map: problem -> list_of_schemes
global_scheme_map = {		
					"TC" 	: ["BF","1H","REDEX","PGP","BSP1","BSP2","DOULION","COLORFUL","BASE"],
					"JP-JC" : ["BF","1H","REDEX","PGP","BASE"],
					"JP-CN" : ["BF","1H","BASE"],
					"JP-OV" : ["BF","1H","BASE"],
					"4C" 	: ["BF","1H","BASE"],
}

# Set threshold parameters we want show for a given scheme
# TODO: Make sure we use the right ones here
global_threshold_map = {	"BF" 		: 0.5,
							"1H" 		: 0.01,
							"REDEX" 	: 0.5,		#0.2 0.5 0.8
							"PGP" 		: 0.5,		#0.2 0.5 0.8
							"BASE" 		: 0,
							"DOULION" 	: 0,
							"COLORFUL" 	: 0,
							"BSP1"		: 50,		#2 5 10
							"BSP2"		: 5,		#2 5 10
			 }


# Adapt names of some graphs
def adapt_graph_names(graph_names):
	graph_names_fixed = []
	name_map = {
		"Si10H16.el" 		: "ch-Si10H16.el",
		"SiO.el" 			: "ch-SiO.el",
		"ca-AstroPh.el" 	: "int-citAsPh.el",
		"gupta3.el" 		: "sc-OptGupt.el",
		"ted_AB.el" 		: "sc-ThermAB.el",
		"p-hat1500-3.el" 	: "dimacs-hat1500-3.el",
		}
	for name in graph_names:
		if name in name_map:
			graph_names_fixed.append(name_map[name])
		else:
			graph_names_fixed.append(name)
	return graph_names_fixed

def read_data(filename, vertex_cnt_filter = None, graphs = None):
	print("Reading file " + filename )
	df = pd.read_csv(filename);
	if graphs != None:
		df = df[(df["graph-name"].isin(graphs))]
	if vertex_cnt_filter != None:
		df = df[(df["vertices"]==vertex_cnt_filter)]

	graph_list = list(df.sort_values(["vertices"], ascending = False)["graph-name"].unique())

	# Collect exact values for each graph and use them to create the ratio column
	true_counts = {};
	base_time = {};

	for g in graph_list:
		try:
			if vertex_cnt_filter == None:
				true_counts[g] = int((df[(df["graph-name"]==g) & (df["approximation-scheme"]=="BASE")])['approximated-count'].values[0]);
				base_time[g] = float((df[(df["graph-name"]==g) & (df["approximation-scheme"]=="BASE")])['total-runtime'].values[0]);
			else:
				true_counts[g] = int((df[(df["graph-name"]==g) & (df["approximation-scheme"]=="BASE")])['approximated-count'].values[0]);
				base_time[g] = float((df[(df["graph-name"]==g) & (df["approximation-scheme"]=="BASE")])['total-runtime'].values[0]);
		except:
			print("Unable to extract exact count and baseline runtime")
			base_time[g] = -1;
			true_counts[g] = -1;

	df["count_ratio"] = df.apply(lambda x: float(x['approximated-count'])/float(true_counts[x['graph-name']]), axis=1)
	df["mem_ratio"] = df.apply(lambda x: float(x['total-size'])/float(x['CSR-size']), axis=1)
	df["time_vs_baseline"] = df.apply(lambda x: x['total-runtime']*1./base_time[x['graph-name']], axis=1)
	df["time_vs_baseline_count_only"] = df.apply(lambda x: x['tc-time']*1./base_time[x['graph-name']], axis=1)
	df["time_vs_baseline_prep_only"] = df.apply(lambda x: x['preprocessing-time']*1./base_time[x['graph-name']], axis=1)
	df["speed_up"] = df.apply(lambda x: base_time[x['graph-name']]*1./x['total-runtime'], axis=1)
	df["speed_up_count_only"] = df.apply(lambda x: base_time[x['graph-name']]*1./x['tc-time'], axis=1)

	return df



def create_count_vs_speedup_plot(problem, df, plotdir):
	# Set style parameters
	plt.rcParams['text.usetex'] = global_latex_en
	plt.rcParams['font.size'] = 24

	# Determine what to plot + labels
	x = 'speed_up_count_only'
	xlabel = "Speed-up vs baseline"
	y = 'count_ratio'
	ylabel = "Relative count"

	# Determine output file
	savename = "count_vs_speedup_" + problem + ".pdf"

	# Create plot
	fig = plt.figure()
	fig, ax = plt.subplots(1,1, figsize = (6,6))
	plt.subplots_adjust(left = 0.1, top = 0.8, bottom = 0.15, right = 0.85)
	
	# Grid
	plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha = 0.4)
	plt.minorticks_on()
	plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)


	# Iterate through schemes
	ymax = 0
	for scheme in global_scheme_map[problem]:
		threshold = global_threshold_map[scheme]
		# Query data
		q = "Problem==@problem and `approximation-scheme`==@scheme and " + x + ">0 and " + y +" >0 and treshold==@threshold"
		tmp_df = df.query(q).sort_values(by=x)
		xlist  = np.array(tmp_df[x]); 
		ylist  = np.array(tmp_df[y]);
		# Data not found
		if len(xlist) == 0 or len(ylist) == 0:
			print("Data not found: %s" + str((problem, scheme, threshold)))
			continue
		# Settings
		ymax = max(ymax, np.max(ylist))
		clist = np.array(tmp_df["mem_ratio"])
		s = 200 if scheme == "BASE" else 100
		# Plot
		sc = ax.scatter(	x = xlist, y = ylist, c = clist, s=s, 
							label = global_label_map[scheme], 
							marker = global_marker_map[scheme], 
							linewidth=1,
							cmap = "Greys", edgecolor = "black", 
							vmin = 1, vmax = 1.5, alpha = 0.8)

	# Set axis size
	if ymax < 1.5: ymax = 1.5;
	if ymax > 5: ymax = 5;
	ax.set_ylim(bottom = 0, top = ymax);
	#ax.set_aspect(1./ax.get_data_ratio())

	# Set axis labels
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	
	# Add color bar legend
	plt.colorbar(sc).set_label(label="Relative memory")

	# Add label legend
	plt.legend(	prop={"size":18} ,ncol = 1, 
				labelspacing = 0.1, handletextpad=0.02, 
				loc = "upper right", bbox_to_anchor=(1, 1), 
				borderpad = 0.1, columnspacing = 0.1, handlelength = 1.5)
	
	# Store plot
	plt.savefig(plotdir + savename, format = 'pdf', dpi=300, bbox_inches = "tight")
	plt.close();


def create_bar_plot(problem, graphs, df, plotdir, plt_width):
	# Set style parameters
	plt.rcParams['text.usetex'] = global_latex_en
	plt.rcParams['font.size'] = 24

	# Determine output file
	savename = "bar_plot_" + problem + ".pdf"

	# Prepare three subplots
	measures = ["speed_up_count_only", "count_ratio", "mem_ratio"];
	names = {mes : n for mes,n in zip(measures, ["   SPEED-UP   ","RELATIVE COUNT","RELATIVE MEMORY"])}
	values = {mes : {} for mes in measures};

	# Prepare data structure to stare data later
	schemes = [s for s in global_scheme_map[problem]]
	for mes in measures:
		for scheme in schemes:
			values[mes][scheme] = {g : -1 for g in graphs};

	used_graphs = []

	# Read data
	q = "Problem==@problem";
	for graph in graphs:				
		graph_ok = True;
		tmp_q = q + " and `graph-name`==@graph"

		for mes in measures:
			values[mes]["BASE"][graph] = 1.

		for scheme in schemes:
			if (scheme != "BASE"):
				if (scheme != "COLORFUL" and scheme != "DOULION"):
					threshold = global_threshold_map[scheme];		 
					tmp_q_2 = tmp_q + " and `approximation-scheme`==@scheme and `treshold`==@threshold";
				else:
					tmp_q_2 = tmp_q + " and `approximation-scheme`==@scheme";
				try:
					for mes in measures:
						val = df.query(tmp_q_2)[mes].values[0];
						values[mes][scheme][graph] = val;
						if val < 0:
							raise ValueError('Negative value found for mes' + mes + " and scheme " + scheme)
				except Exception as e:
					print("Exception while reading data: %s %d %s" % (scheme, threshold, graph))
					print(e)
					graph_ok = False
		if graph_ok:		
			used_graphs.append(graph)

	# Create plot
	plt.style.use('grayscale')	
	fig, axs = plt.subplots(3, figsize = (plt_width,7), sharex= True)
	plt.subplots_adjust(wspace=0, hspace=0.2)

	# Configure axis
	colormap = plt.cm.Greys
	for ax in axs:
		ax.tick_params(axis='x', colors='grey')
		ax.grid(b=True, which='major', color='#999999', linestyle='-', alpha=0.3, linewidth=0.6)
		ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.4)
		ax.set_prop_cycle(color=[colormap(i) for i in np.linspace(0, 1,len(schemes))])

	# Adapt graph names to what we want to show
	graph_names_adapted = adapt_graph_names(used_graphs)
	graph_names = [g.split(".")[0].replace("_","\_") for g in graph_names_adapted]

	# Format x-axis 
	positions = np.array(range(len(used_graphs)));
	space = 0.2;
	width = (1. - space)/len(schemes);


	# Create bar plot
	for ax, mes in zip(axs, measures):
		rel_pos = (-1. + space)/2
		for scheme in schemes:
			label = (global_label_map[scheme] if scheme in global_label_map else scheme)
			ax.bar(positions + rel_pos,[values[mes][scheme][g] for g in used_graphs], width = width, label = label, edgecolor = "black", lw = 0.3);
			ax.set_title(names[mes], fontsize = 18, ha = "left", x = 0, y = 0.975)
			rel_pos += width;
	

	# x-ticks and graph names
	for ax in axs[:-1]:
		ax.xaxis.set_visible(False)
	axs[-1].set_xticks(positions);
	axs[-1].minorticks_off();
	axs[-1].set_xticklabels(graph_names, color = "black");
	axs[-1].xaxis.set_visible(True)

	# y-ticks
	(ymin, ymax) = axs[1].get_ylim()
	for ax in axs:
		(ymin, ymax) = ax.get_ylim()
		positions = [x * math.ceil(ymax) / 4 for x in range(5)]
		ax.set_yticks(positions);
		ax.set_xlim(-0.6,20.5)
	fig.autofmt_xdate(rotation=20)  
	axs[1].set_ylim(0,1.3)

	# legend
	h, l = axs[0].get_legend_handles_labels() # Extracting handles and labels
	ph = [plt.plot([],marker="", ls="")[0]] # Canvas
	plt.legend(h, l, loc='upper left', ncol = int(plt_width/2), frameon=False,prop={"size":22.5},bbox_to_anchor=(0.0, 3.45), columnspacing = 0.5, handlelength = 0.5)

	# Save result
	plt.savefig(plotdir + savename, format = 'pdf', dpi=300, bbox_inches = "tight")
	plt.close();

def create_all_plots():

	graph_list	= [
		"SiO.el",
		"ca-AstroPh.el",
		"Si10H16.el",
#		"gupta3.el",
		"bio-WormNet-v3.el",
		"bio-CE-CX.el",
		"ted_AB.el",
		"bio-HS-CX.el",
		"bio-HS-LC.el",
		"bio-DM-CX.el",
		"bio-DR-CX.el",
		"econ-psmigr1.el",
		"econ-psmigr2.el",
		"econ-orani678.el",
		"bio-SC-HT.el",
		"bio-CE-PG.el",
		"bio-SC-GT.el",
		"p-hat1500-3.el",
		"econ-beaflw.el",
		"econ-beacxc.el",
		"econ-mbeacxc.el",
		"bn-mouse_brain_1.el",
	]

	problem_list = ["TC","JP-JC","JP-CN","JP-OV"]

	graph_list_debug = ["kronecker"]

	#vertex_cnt_filter = 1048576
	vertex_cnt_filter = None 


	for problem in problem_list:
		filename = "../../real_graph_results/" + problem.lower() + "-real.csv"
		df = read_data(filename, vertex_cnt_filter, graph_list)
		#create_count_vs_speedup_plot(problem, df, "plots/")
		plt_width = 26 if problem == "TC" else 12
		create_bar_plot(problem, graph_list, df, "real_graphs_images/", plt_width)


create_all_plots()


import matplotlib
import csv
import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mtick
import os
import sys
import statistics

global_colors = ['#CCCCCC', '#999999', '#666666', '#333333', '#000000']
global_estimators = {
	"single_set" : [	r'$\widehat{|X|}_L$', 
						r'$\widehat{|X|}_S$', 
						r'$\widehat{|X|}_?$',	
					],
	"intersection": [	r'$|\widehat{X \cap Y}|_{AND}$', 
						r'$|\widehat{X \cap Y}|_{OR}$', 
						r'$|\widehat{X \cap Y}|_{L}$',
						r'$|\widehat{X \cap Y}|_{1H}$',
						r'$|\widehat{X \cap Y}|_{kH}$',
					],
					}

def adapt_graph_names(graph_names):
	graph_names_fixed = []
	name_map = {
		"Si10H16" : "ch-Si10H16",
		"SiO" : "ch-SiO",
		"ca-AstroPh" : "int-citAsPh",
		"gupta3" : "sc-OptGupt",
		"ted_AB" : "sc-ThermAB",
		"p-hat1500-3" : "dimacs-hat1500-3",
		}
	for name in graph_names:
		if name in name_map:
			graph_names_fixed.append(name_map[name])
		else:
			graph_names_fixed.append(name)
	return graph_names_fixed

def round_up(n, decimals=0):
	multiplier = 10 ** decimals
	return math.ceil(n * multiplier) / multiplier

def read_config(config_file):
	with open("./intersection_estimator_config/" + config_file + ".csv", 'r') as csvfile:
		csvreader = csv.reader(csvfile)
		graph_names = []
		for row in csvreader:
			graph_names.append(row[0])
		return graph_names

def read_graph(graph_file, skip = False):
	with open("./intersection_estimator_results/" + graph_file + ".csv", 'r') as csvfile:
		csvreader = csv.reader(csvfile)
		estimators = []
		for row in csvreader:
			correct = int(row[0])			
			if correct == 0 or min([int(row[i]) for i in range(len(row))]) < 0:
				continue
			for i in range(len(row) - 1):
				if len(estimators) <= i:
					estimators.append([])
				est = int(row[i+1])
				metric = abs(est-correct)/correct
				estimators[i].append(metric)
		return estimators

def create_plots_with_fixed_parameters(config_file, b_values, mem_values, estimator_type, metric, ylim):
	graph_names = read_config(config_file)
	for b in b_values:
		for mem in mem_values:
			create_plot_with_fixed_parameters(graph_names, estimator_type, b, mem, metric, ylim)

def create_plots_with_varying_b(config_file, b_values, mem_values, estimator_type, metric):
	graph_names = read_config(config_file)
	for mem in mem_values:
		create_plot_with_varying_b(graph_names, estimator_type, b_values, mem, metric)

def create_plots_with_varying_mem(config_file, b_values, mem_values, estimator_type, metric):
	graph_names = read_config(config_file)
	for b in b_values:
		create_plot_with_varying_mem(graph_names, estimator_type, b, mem_values, metric)

def create_plot_with_fixed_parameters(graph_names, estimator_type, b, mem, metric, ylim):
	n_graphs = len(graph_names)
	data = []
	# Read all data necessary
	used_graph_names = []
	for i in range(n_graphs):
		graph_name = graph_names[i]
		graph_file = estimator_type + "_" + graph_name + "_" + str(b) + "_" + str(mem)
		estimators = read_graph(graph_file)	
		if len(estimators) == 0:
			continue			
		n_est = len(estimators)
		for estimator in estimators:
			data.append(estimator)
		data.append([])
		used_graph_names.append(graph_name)
	graph_names = used_graph_names
	n_graphs = len(graph_names)
	
	# Initialize plot
	plt.close("all")
	plt.rcParams['text.usetex'] = True
	width = 12 if n_graphs > 8 else 6
	fig, ax = plt.subplots(1,1, figsize = (width,8))
	left = 0.13 if estimator_type == "single_set" else 0.25
	bottom = 0.25 if n_graphs > 3 else 0.1
	plt.subplots_adjust(left=left, right = 0.99, bottom = bottom, top = 0.95)
	# Create plot
	xvalues = [i for i in range((n_est + 1) * n_graphs) if i % (n_est + 1) != n_est]
	bp = ax.boxplot(data, xvalues, showfliers=False, patch_artist=True, widths = 0.8)
	colors = (global_colors[:n_est] + ["blue"]) * n_graphs
	for patch, color in zip(bp['boxes'], colors):
		patch.set_facecolor(color)
	# Format x-axis
	graph_names_fixed = adapt_graph_names(graph_names)
	xlabels = [graph_names_fixed[i].replace("_","\_") for i in range(n_graphs)]
	center = (n_est + 1) / 2
	xlocations = [(n_est + 1) * i + center for i in range(n_graphs)]
	rotation = 45 if n_graphs > 3 else 0
	ha = "right" if n_graphs > 3 else "center"
	plt.xticks(xlocations, xlabels, rotation=rotation, ha = ha, fontsize = 24, fontname="sans-serif")
	ax.set_xlim(-1, (n_est + 1) * n_graphs + 1)
	# Format y-axis
	ax.set_ylim(-0.05,ylim)
	ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, symbol='%', is_latex=False))
	plt.yticks(fontsize = 24)
	ytext = "Relative difference: " + r"$\frac{|\widehat{"+metric+"}-"+metric+"|}{"+metric+"}$"
	ax.set_ylabel(ytext, fontsize = 24)
    # Create legend 
	custom_lines = [Line2D([0], [0], color=colors[i], lw=4) for i in range(n_est)]
	estimators = global_estimators[estimator_type][:n_est]
	leg = ax.legend(custom_lines, estimators, loc = "upper center", fontsize = 20, ncol = 2, columnspacing = 1, handletextpad = 0.2)
	# Misc
	ax.grid(axis = "y")
	# Save figure as pdf
	plt.savefig("./intersection_estimator_plots/plot_" + estimator_type + "_b" + str(b) + "_mem" + str(mem) + ".pdf")

def create_plot_with_varying_b(graph_names, estimator_type, b_values, mem, metric):
	n_graphs = len(graph_names)
	n_bs = len(b_values)
	data = []
	# Read all data necessary
	for i in range(n_graphs):
		graph_name = graph_names[i]
		for b in b_values:
			graph_file = estimator_type + "_vb_" + graph_name + "_" + str(b) + "_" + str(mem)
			estimators = read_graph(graph_file)	
			for estimator in estimators:
				data.append(estimator)
			data.append([])
	n_est = len(estimators)
	# Initialize plot
	plt.close("all")
	plt.rcParams['text.usetex'] = True
	fig, ax = plt.subplots(1,1, figsize = (6,6))
	plt.subplots_adjust(left=0.25, right = 0.99, bottom = 0.13, top = 0.97)
	# Create plot
	xvalues = [i for i in range((n_est + 1) * n_bs * n_graphs) if i % (n_est + 1) != n_est]
	bp = ax.boxplot(data, xvalues, showfliers=False, patch_artist=True, widths = 0.8)
	colors = (global_colors[:n_est] + ["blue"]) * n_graphs * n_bs
	for patch, color in zip(bp['boxes'], colors):
		patch.set_facecolor(color)
	# Format x-axis
	xlabels = []
	xlocations = []
	for i in range(n_graphs):
		for j in range(n_bs):
			xlabels.append(str(b_values[j]))	
			xlocations.append(i * n_bs * (n_est + 1) + j * (n_est + 1) + (n_est + 1) / 2)
	plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=True)
	plt.xticks(xlocations, xlabels, fontsize = 24)
	ax.set_xlabel("Hash function count b", fontsize = 24)
	ax.set_xlim(-1, (n_est + 1) * n_graphs * n_bs + 1)
	# Format y-axis
	ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, symbol='%', is_latex=False))
	plt.yticks(fontsize = 24)
	ytext = "Relative difference: " + r"$\frac{|\widehat{"+metric+"}-"+metric+"|}{"+metric+"}$"
	ax.set_ylabel(ytext, fontsize = 24)
    # Create legend 
	custom_lines = [Line2D([0], [0], color=colors[i], lw=4) for i in range(n_est)]
	estimators = global_estimators[estimator_type][:n_est]
	leg = ax.legend(custom_lines, estimators, loc = "upper left", fontsize = 24)
	# Misc
	ax.grid(axis = "y")
	# Save figure as pdf
	plt.savefig("./intersection_estimator_plots/plot_" + estimator_type + "_vary_b" + "_mem" + str(mem) + ".pdf")

def create_plot_with_varying_mem(graph_names, estimator_type, b, mem_values, metric):
	n_graphs = len(graph_names)
	n_mems = len(mem_values)
	data = []
	# Read all data necessary
	for i in range(n_graphs):
		graph_name = graph_names[i]
		for mem in mem_values:
			graph_file = estimator_type + "_vm_" + graph_name + "_" + str(b) + "_" + str(mem)
			estimators = read_graph(graph_file)	
			for estimator in estimators:
				data.append(estimator)
			data.append([])
	n_est = len(estimators)
	# Initialize plot
	plt.close("all")
	plt.rcParams['text.usetex'] = True
	fig, ax = plt.subplots(1,1, figsize = (6,6))
	plt.subplots_adjust(left=0.25, right = 0.99, bottom = 0.13, top = 0.97)
	# Create plot
	xvalues = [i for i in range((n_est + 1) * n_mems * n_graphs) if i % (n_est + 1) != n_est]
	bp = ax.boxplot(data, xvalues, showfliers=False, patch_artist=True, widths = 0.8)
	colors = (global_colors[:n_est] + ["blue"]) * n_graphs * n_mems
	for patch, color in zip(bp['boxes'], colors):
		patch.set_facecolor(color)
	# Format x-axis
	xlabels = []
	xlocations = []
	for i in range(n_graphs):
		for j in range(n_mems):
			xlabels.append(str(mem_values[j]) + "\%")	
			xlocations.append(i * n_mems * (n_est + 1) + j * (n_est + 1) + (n_est + 1) / 2)
	plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=True)
	plt.xticks(xlocations, xlabels, fontsize = 24)
	ax.set_xlabel("Storage budget s", fontsize = 24)
	ax.set_xlim(-1, (n_est + 1) * n_graphs * n_mems + 1)
	# Format y-axis
	ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, symbol='%', is_latex=False))
	plt.yticks(fontsize = 24)
	ytext = "Relative difference: " + r"$\frac{|\widehat{"+metric+"}-"+metric+"|}{"+metric+"}$"
	ax.set_ylabel(ytext, fontsize = 24)
    # Create legend 
	custom_lines = [Line2D([0], [0], color=colors[i], lw=4) for i in range(n_est)]
	estimators = global_estimators[estimator_type][:n_est]
	leg = ax.legend(custom_lines, estimators, loc = "upper right", fontsize = 24)
	# Misc
	ax.grid(axis = "y")
	# Save figure as pdf
	plt.savefig("./intersection_estimator_plots/plot_" + estimator_type + "vary_mem" + "_b" + str(b) + ".pdf")

create_plots_with_fixed_parameters("graphs_intersection", [1,4], [33], "intersection", "|X \cap Y|", 2.5)


import matplotlib.pyplot as plt
import itertools
import csv
import sys

global_color_map = {	"BASE" : 		"#000000",
						"1H" : 			"#000099",
						"BF" : 			"#009900",
						"COLORFUL" : 	"#990000",
						"DOULION" : 	"#990099",
				}

global_marker_map = {	"BASE" : 		"o",
						"1H" : 			"x",
						"BF" : 			"^",
						"COLORFUL" : 	"*",
						"DOULION" : 	"d",
				}

global_label_map = {
					"BASE" : "Exact TC",
					"BF" : "ProbGraph (BF)",
					"1H" : "ProbGraph (1H)",	
					"DOULION" : "Doulion",		
					"COLORFUL" : "Colorful",		
				}


def read_data(filename, problem):
	estimators = []
	threadss = {}
	runtimes = {}
	edge_densities = {}
	# Open file
	with open("../../scaling_experiments/" + filename + ".csv", newline='') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='"')
		# Get column headers
		columns = list(next(reader))
		# Iterate through rows
		for row in reader:	
			if row[columns.index("Problem")] != problem.upper():
				continue
			# Read data
			estimator = row[columns.index("approximation-scheme")]
			threads = int(row[columns.index("thread-count")])
			runtime = float(row[columns.index("total-runtime")])
			vertices = int(row[columns.index("vertices")])
			edges = int(row[columns.index("edges")])
			edge_density = float(edges) / float(vertices) 
			# Store data 
			if estimator not in estimators:
				estimators.append(estimator)	
				threadss[estimator] = []
				runtimes[estimator] = []
				edge_densities[estimator] = []
			if threads not in threadss[estimator]:
				threadss[estimator].append(threads)
				runtimes[estimator].append(runtime)
				edge_densities[estimator].append(edge_density)
			elif runtime > runtimes[estimator][threadss[estimator].index(threads)]:
				runtimes[estimator][threadss[estimator].index(threads)] = runtime
				edge_densities[estimator].append(edge_density)
	for estimator in estimators:
		edge_densities[estimator] = [x for _,x in sorted(zip(threadss[estimator],edge_densities[estimator]))]
		runtimes[estimator] = [x for _,x in sorted(zip(threadss[estimator],runtimes[estimator]))]
		threadss[estimator] = sorted(threadss[estimator])
	return (estimators, threadss, runtimes, edge_densities)

def create_one_plot(problem, filename, plotname, show_lin_speedup, show_edge_density):
	# Initialize Plot
	fig, ax = plt.subplots(1,1, figsize = (4,4))
	plt.subplots_adjust(left=0.16, right = 0.99, top = 0.98, bottom = 0.125)
	plt.rcParams['text.usetex'] = True
	# Read data
	(estimators, threadss, runtimes, edge_densities) = read_data(filename, problem)
	if len(estimators) == 0:
		return
	# Create plot
	tmp = set()
	for i in range(len(estimators)):
		estimator = estimators[i]
		label = global_label_map[estimator]
		marker = global_marker_map[estimator]
		color = global_color_map[estimator]
		# Plot data series
		if show_lin_speedup:
			ideal = [runtimes[estimator][0] / threadss[estimator][i] for i in range(len(threadss[estimator]))]
			ideal_color = color.replace("99","FF").replace("00","99")
			ax.plot(threadss[estimator], ideal, label = "Lin. Speedup " + label, marker = marker, fillstyle="none", color=ideal_color)
		if show_edge_density:
			ax.text(0.8,1.39, "$\\frac{m}{n}\\!\\!\\approx$" ,ha='center', va='center', fontsize = 13)
			for i in range(len(threadss[estimator])):
				if threadss[estimator][i] not in tmp:
					#txt = "$\\frac{m}{n}\\approx" + str(int(round(edge_densities[estimator][i],0))) +"$"
					txt = str(int(round(edge_densities[estimator][i],0)))
					ax.text(threadss[estimator][i],1.39, txt,ha='center', va='center', fontsize = 13)
					tmp.add(threadss[estimator][i])
			# Hack to fix left and right spacing to make edge density fit
			ax.plot(38,2)
			ax.plot(0.8,2)
		ax.plot(threadss[estimator], runtimes[estimator], label = "" + label, marker = marker, fillstyle='none', color=color)
	# Axis settings
	ax.grid()
	ax.set_ylim(bottom = 1)
	ax.set_xlabel("Number of Threads",fontsize = 14)
	ax.set_ylabel("Runtime [s]", fontsize = 14)
	ax.set_xscale('log', base = 2)
	ax.set_yscale('log', base = 10)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	# Legend
	ncol = 1
	pos = ("lower left",(0,0)) if "strong" in plotname else ("upper left",(0,1))

	handles, labels = plt.gca().get_legend_handles_labels()
	times = [runtimes[estimators[i]][1] for i in range(len(estimators))]
	order = [x for _,x in sorted(zip(times,list(range(len(times)))), reverse = True)]
	plt.legend( [handles[idx] for idx in order],[labels[idx] for idx in order])
	plt.legend(	[handles[idx] for idx in order],
				[labels[idx] for idx in order],
				ncol = ncol,							# Columns and order of entries
				loc=pos[0], bbox_to_anchor=pos[1], 		# Position of legend
				prop={'size': 12.5}, markerscale = 0.9,					# Text and marker size
				handletextpad=0.1, handlelength = 0.9, 					
				columnspacing = 1, labelspacing = 0.2, 
				)
	plt.savefig("../../scaling_images/" + plotname + ".pdf")
	plt.close()

def create_all_plots():
	# Strong scaling
	problems = ["jp-cn","jp-jc","jp-ov","tc"]
	for problem in problems:
		create_one_plot(problem, "strong_scaling_" + problem, "plot_strong_scaling_" + problem, False, False)
	# Weak scaling
	for problem in problems:
		create_one_plot(problem, "weak_scaling_" + problem, "plot_weak_scaling_" + problem, False, True)


create_all_plots()

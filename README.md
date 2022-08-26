### Prerequirisites: #

  * boost 1.76
  * g++ 8.3.0
  * python 3.7

### Building: #

Run the following commands from the main probgraph-public folder:
```
cd src
make
cd ../ 
```

### Using ProbGraph: #

Now run the test to check that everything works:

`./test.sh`

or you can try it yourself by running any of the executables in probgraph-public/src/

the general formulation is

`./src/PROBLEM_APPROXIMATOR OPTIONS`

possible choices for PROBLEM are:
  * `tc` : Triangle Counting
  * `jp-jc` : Clustering (Jaccard)
  * `jp-ov` : Clustering (Overlap)
  * `jp-cn` : Clustering (Common neigh.)
  * `4c`: 4-cliques

possible choices for APPROXIMATOR are:
  * `base` : Baseline
  * `1h`: one-hash
  * `bf`: Bloom Filter
  * `pgp`: Partial Graph Processing
  * `redex`: Reduced Execution
  * `colorful`, `doulion` : only available for Triangle Counting

For example, 

`./src/tc_bf -g 17 -k 256`

will run Triangle counting (TC) using Bloom Filters (BF) on a Kronecker graph with 2^17 nodes and 256 avg. degree (-g 17 -k 256).
Run any executable with the -h flag to see the available options

In order to run experiments on stored graphs, you will have to download them in the probgraph-sc/src/graphs/ folder. Some small graphs are already in the folder. You can find a list of graphs in the paper together with references to the datasets. In alternative, you can download your own graphs in the graphs folder. Be sure they are in the ordered edgelist format. 

### Reproducing our results: #

To reproduce the experiments in the paper, run 

`./launch_experiments` (use the _CLUSTER version if you are on a cluster) 

Then, collect all results using 

`./collect_results.sh`

this will overwrite the existing csv files in the result folders. If you just want to generate the images from the existing data, skip the two previous steps and just run

`./generate_images.sh`

### Sample images: #

 Analysis of performance, accuracy and memory for the BF and MH ProbGraph approximators (Jaccard Similarity Clustering) 
 ![Analysis of performance, accuracy, and memory of ProbGraph](/sample_images/barplot_test_real_JP-JC.pdf)

Accuracy and Speed-up for various measures estimated with ProbGraph on Kronecker graphs.
![Advantages of ProbGraph for Kronecker graphs.](/sample_images/main-results-kron___low-mem.pdf)

Weak scaling behaviour of ProbGraph and some representative baselines for the TC problem.
![Weak scaling for the TC problem](/sample_images/plot_weak_scaling_tc.pdf)


Analysis of the accuracy of ProbGraph estimators for intersection.
![Analysis of the accuracy of ProbGraph estimators for intersection.](/sample_images/plot_intersection_b4_mem33.pdf)

### Reproducing the results on intersection estimator accuracy: #

Run the following commands to reproduce our accuracy study for the interesection estimator: 

`cd src/src`

`g++ -o evaluate_intersection_estimators evaluate_intersection_estimators.cpp MurmurHash3.cpp`

`./evaluate_intersection_estimators`

`cd ../..`

`python3 create_intersection_estimator_plots.py`

The plots will be saved in the folder 'intersection_estimator_plots'. Please note that due to the random choice of hash functions, the plots will slightly vary each time they are re-generated.

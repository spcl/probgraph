// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <memory>
#include "../benchmark.h"
#include "../builder.h"
#include "../command_line.h"
#include <sys/time.h>
#include "../graph.h"
#include "../pvector.h"
#include "../set_operation.hpp"
#include "../sets.hpp"
/*
GAP Benchmark Suite
Kernel: Triangle Counting (TC)
Author: Scott Beamer

Will count the number of triangles (cliques of size 3)

Requires input graph:
  - to be undirected
  - no duplicate edges (or else will be counted as multiple triangles)
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

This implementation reduces the search space by counting each triangle only
once. A naive implementation will count the same triangle six times because
each of the three vertices (u, v, w) will count it in both ways. To count
a triangle only once, this implementation only counts a triangle if u > v > w.
Once the remaining unexamined neighbors identifiers get too big, it can break
out of the loop, but this requires that the neighbors to be sorted.

Another optimization this implementation has is to relabel the vertices by
degree. This is beneficial if the average degree is high enough and if the
degree distribution is sufficiently non-uniform. To decide whether or not
to relabel the graph, we use the heuristic in WorthRelabelling.
*/


using namespace std;
inline NodeID* get_vcandidates_start(const Graph &g, const NodeID u) {
  NodeID* uneigh = g.out_neigh(u).begin();
  NodeID nuneigh = g.out_degree(u);
  NodeID* vcandidates_start;
  if (nuneigh==0) {
    vcandidates_start = uneigh;
  } else {
    int oldind = 0;
    int ind = 0;
    while (ind < nuneigh) {
      if (*(uneigh+ind)>u) {
        break;
      } else {
        oldind = ind;
        ind = ((ind+1)<<1)-1; // 2x(ind+1) - 1
      }
    }
    long int* low = uneigh + oldind;
    long int* high;
    if (ind<nuneigh) {
      high = uneigh + ind;
    }
    else {
      if (u > *(uneigh + nuneigh-1)) {
        high = uneigh + nuneigh;
      } else {
        high = uneigh + nuneigh-1;
      }
    }
    while (high>low) {
      long int* pivot = (low + (high-low)/2);
      if (*pivot<u) low = pivot + 1;
      else high = pivot;
    }
    vcandidates_start = high;
  }
  return vcandidates_start;
}

size_t OrderedCount(const Graph &g, float threshold_, int k_, FILE *fp_tr, std::string graphName, int threads) {
  double tc_time = -1;
  double pp_time = -1;
  double approx_str_size = -1;
  double initial_csr_size = -1;
  int num_thresholds;
  float *thresholds;
  if (fp_tr != nullptr) {
    fscanf(fp_tr, "%d", &num_thresholds);
    printf("Running for %d threshold values\n", num_thresholds);
    assert (num_thresholds >= 1);
    thresholds = (float *) malloc (sizeof(float) * num_thresholds);
    for (int i = 0; i < num_thresholds; i++) 
      fscanf(fp_tr, "%f", &thresholds[i]);
  }
  else {
    num_thresholds = 1;
    thresholds = (float *) malloc (sizeof(float) * num_thresholds);
    thresholds[0] = threshold_;
  }
  
  for (int i = 0; i < num_thresholds; i++) {
    double threshold = thresholds[i];
    printf("-------------------------------\n");
    printf("Running for threshold: %f\n", threshold);
    Timer t;
    t.Start();
    struct timeval time;
    gettimeofday(&time, NULL);
    int seed = time.tv_usec;
    // One OneHashSet per node (neighborhood)
    vector<OneHashSet<NodeID>*> sets;
    sets.reserve(g.num_nodes());
    vector<NodeID*> hstart;
    hstart.reserve(g.num_nodes());
    size_t total = 0;
    int64_t max_degree = 0;

    // Find max degree and use as value for k in combination with threshold
    #pragma omp parallel for reduction(max:max_degree) schedule(static, 32)
    for(NodeID u=0; u<g.num_nodes(); ++u) {
      max_degree = std::max(g.out_degree(u), max_degree);
    }

    int k = k_;
    if (k <= 0) k = max_degree * threshold;
    if (k < 1) k = 1;
    std::cout << "max deg: " << max_degree << std::endl;
    std::cout << "k: " << k << std::endl;

    #pragma omp parallel for schedule(dynamic, 64)
    for(NodeID u = 0; u<g.num_nodes(); ++u) {
      hstart[u] = get_vcandidates_start(g,u);
      OneHashSet<NodeID> * neigh = new OneHashSet<NodeID>(g.out_neigh(u).begin(), hstart[u], k, seed);
      sets[u] = neigh;
    }
    t.Stop();
    PrintTime("Datastructure building", t.Seconds());
    pp_time = t.Seconds();

    size_t size = 0;
    for(int i = 0; i<g.num_nodes(); i++){
      size+= sets[i]->total_size();
    }
    PrintTime("Datastructure size", size);

    approx_str_size = size / (1024.0 * 1024.0); // MB
    initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

    t.Start();
    #pragma omp parallel reduction(+:total)
    {
      #pragma omp for schedule(dynamic, 64)
      for (NodeID u = 0; u < g.num_nodes(); ++u) {
        for (NodeID* vp = g.out_neigh(u).begin(); vp < hstart[u]; vp++) {
          int n = (*sets[u]).intersect_count(*sets[*vp]);
          total += n;
        }
      }
    }
    t.Stop();
    //  PrintTime("TC runtime", t.Seconds());
    PrintTime("Intersection time", t.Seconds());
    tc_time = t.Seconds();

    #pragma omp for schedule(dynamic, 64)
    for(size_t i = 0; i<sets.size(); i++){
      delete sets[i];
    }
    sets.clear();
    sets.shrink_to_fit();
    hstart.clear();
    hstart.shrink_to_fit();
    std::cout << "triangles: " << total << std::endl;
    printf("-------------------------------\n");

    int m = 0;

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (1-Hash parameter)] [k (another 1-Hash parameter-number of min hashes] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR TC 1H TC_1H " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << pp_time << " " << tc_time << " " << pp_time + tc_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (1-Hash parameter)] [k (another 1-Hash parameter-number of min hashes] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS TC 1H TC_1H " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;

    printf("ooo -------------------------------\n");
  }
  free(thresholds);
  // return total;
  return 0;
}


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
  cout << total_triangles << " triangles" << endl;
}

// Compares with simple serial implementation that uses std::set_intersection
bool TCVerifier(const Graph &g, size_t test_total, const CLApp &cli) {
  size_t total = 0;
  #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    for (NodeID v : g.out_neigh(u)) {
      if (v > u)
        break;
      auto it = g.out_neigh(u).begin();
      for (NodeID w : g.out_neigh(v)) {
        if (w > v)
          break;
        while (*it < w)
          it++;
        if (w == *it)
          total++;
      }
    }
  }

  cout << "approx total: " << test_total << endl;
  cout << "exact total: " << total << endl;
  cout << "acc: " << (float)((float)test_total-total)/total*100  << "\% error: true " << total << " counted " << test_total << endl;
  if (total != test_total)
    cout << total << " != " << test_total << endl;
  return total == test_total;
}

int main(int argc, char* argv[]) {
  CLSIMDApp cli(argc, argv, "kMinHash triangle count");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  auto ct = [&cli] (const Graph& g){
     return OrderedCount(g, cli.treshold(), cli.k(), cli.get_parameters_file(), cli.getGraphBasename(), cli.getThreadNum());
   };
  BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);
  return 0;
}

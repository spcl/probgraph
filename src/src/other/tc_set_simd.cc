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

size_t OrderedCount(const Graph &g, const vector<KMinHashSet<NodeID>*> &sets, const vector<NodeID*> &hstart) {
  size_t total = 0;
  #pragma omp parallel reduction(+:total)
  {
  #pragma omp for schedule(dynamic, 16	) 
    for (NodeID u = 0; u < g.num_nodes(); ++u) {
     for (NodeID* vp = g.out_neigh(u).begin(); vp < hstart[u]; vp++) {
        int n = (*sets[g.num_nodes()+u]).intersect_count(*sets[g.num_nodes() + *vp]);
       total += n;
      }
    }
  }
  
  return total;
}


// heuristic to see if sufficently dense power-law graph
bool WorthRelabelling(const Graph &g) {
  int64_t average_degree = g.num_edges() / g.num_nodes();
  if (average_degree < 10)
    return false;
  SourcePicker<Graph> sp(g);
  int64_t num_samples = min(int64_t(1000), g.num_nodes());
  int64_t sample_total = 0;
  pvector<int64_t> samples(num_samples);
  for (int64_t trial=0; trial < num_samples; trial++) {
    samples[trial] = g.out_degree(sp.PickNext());
    sample_total += samples[trial];
  }
  sort(samples.begin(), samples.end());
  double sample_average = static_cast<double>(sample_total) / num_samples;
  double sample_median = samples[num_samples/2];
  return sample_average / 1.3 > sample_median;
}


// uses heuristic to see if worth relabeling
//size_t Hybrid(const Graph &g) {
//  if (WorthRelabelling(g))
//    return OrderedCount(Builder::RelabelByDegree(g));
//  else
//    return OrderedCount(g);
//}


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
  cout << total_triangles << " triangles" << endl;
}


// Compares with simple serial implementation that uses std::set_intersection
bool TCVerifier(const Graph &g, size_t test_total, const CLApp &cli) {
  size_t total = 0;
  vector<NodeID> intersection;
  intersection.reserve(g.num_nodes());
  for (NodeID u : g.vertices()) {
    for (NodeID v : g.out_neigh(u)) {
      auto new_end = set_intersection(g.out_neigh(u).begin(),
                                      g.out_neigh(u).end(),
                                      g.out_neigh(v).begin(),
                                      g.out_neigh(v).end(),
                                      intersection.begin());
      intersection.resize(new_end - intersection.begin());
      total += intersection.size();
    }
  }
  total = total / 6;  // each triangle was counted 6 times
  if (total != test_total)
    cout << total << " != " << test_total << endl;
  return total == test_total;
}

std::vector<NodeID*> addSets(vector<KMinHashSet<NodeID>*> &sets,const Graph &g,float tres){
  for(NodeID u = 0; u<g.num_nodes(); u++){
      KMinHashSet<NodeID> * neigh = new KMinHashSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end(), tres);

      sets.push_back(neigh);

  }
  std::vector<NodeID*> hstart;
  hstart.reserve(g.num_nodes());
  for (NodeID u = 0; u < g.num_nodes(); ++u) {
    NodeID* uneigh = g.out_neigh(u).begin();
    NodeID nuneigh = g.out_degree(u);    
    NodeID* vcandidates_start;
    if(nuneigh==0){
        vcandidates_start=uneigh;
    }else{
        int oldind = 0;
        int ind = 0; 
        while(ind < nuneigh){
            if (*(uneigh+ind)>u){
                break;
            }else{
                oldind = ind;
                ind = ((ind+1)<<1)-1;
            }
        }
        int* low = uneigh + oldind;
        int* high;
        if(ind<nuneigh)
            high = uneigh + ind;
        else{
            if (u > *(uneigh + nuneigh-1))
                high = uneigh + nuneigh;
            else
                high = uneigh + nuneigh-1;
        }
        while(high>low){
            int* pivot = (low + (high-low)/2);
            if (*pivot<u) low = pivot + 1;
            else high = pivot;
        }
        vcandidates_start = high;
                
    }
    
    KMinHashSet<NodeID> * neigh = new KMinHashSet<NodeID>(g.out_neigh(u).begin(), vcandidates_start, tres);
    
    sets.push_back(neigh);
    hstart.push_back(vcandidates_start);

  }
  return hstart;
 
}

int main(int argc, char* argv[]) {
  CLSIMDApp cli(argc, argv, "triangle count");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  vector<KMinHashSet<NodeID>*> sets;
  auto hstart = addSets(sets, g, cli.treshold());
  auto ct = [&sets, &hstart] (const Graph& g){
     return OrderedCount(g, sets, hstart);
   };
  BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);
  for(auto set : sets) delete set;
  return 0;
}

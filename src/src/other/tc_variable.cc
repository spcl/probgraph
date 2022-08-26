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
#include <sys/time.h>
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

// inline NodeID* get_vcandidates_start(const Graph &g, const NodeID u) {

// }

std::vector<NodeID*> addSets(vector<BloomSet<NodeID>*> &sets, const Graph &g, float tres, int k_, int m_) {
    cout << "threshold: " << tres << endl;
    int64_t max_degree = 0;
    for (NodeID u=0; u<g.num_nodes(); ++u) {
      max_degree = std::max(max_degree, g.out_degree(u));
    }
    cout << "max_degree: " << max_degree << endl;

    if (k_ <= 0 and m_ <= 0 ) std::cout<<"k and m will set automatically"<<std::endl;

    int m = m_, k=k_;

    if (m <= 0){
        m = max_degree * tres;
    }
    if (m < 32) {
      m = 32;
    }

    if (k <= 0){
        if (max_degree!=0) {
            k = m/max_degree*0.69314718056; // m/n*ln(2)
        }
    }
    if (k < 1) {
        k = 1;
    }

    cout << "m: " << m << endl;
    cout << "k: " << k << endl;
    struct timeval seed;
    gettimeofday(&seed, NULL);
    srand(seed.tv_usec);
    std::vector<uint32_t> seeds;
    for (int i = 0; i<k; i++) {
      seeds.push_back(rand());
    }


  //for(NodeID u = 0; u<g.num_nodes(); u++){
  //    BloomSet<NodeID> * neigh = new BloomSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end(), m, k);

  //    sets.push_back(neigh);

  //}
  std::vector<NodeID*> hstart; // Stores the end iterator to the range of neighbor vertices
  hstart.reserve(g.num_nodes());
  for (NodeID u = 0; u < g.num_nodes(); ++u) {
    NodeID* u_neigh = g.out_neigh(u).begin();
    int64_t u_deg = g.out_degree(u);    
    NodeID* vcandidates_start; // The end iterator to the range of vertices added to the current bloomSet
    if (u_deg==0) {
      vcandidates_start=u_neigh; // No vertices to be added to bloomset
    } else {
      int oldind = 0;
      int ind = 0; 
      while (ind < u_deg) {
        if (*(u_neigh+ind)>u) {
          break;
        } else {
          oldind = ind;
          ind = ((ind+1)<<1)-1; // 2x(ind+1) - 1
        }
      }
      int* low = u_neigh + oldind;
      int* high;
      if (ind<u_deg) {
        high = u_neigh + ind;
      }
      else {
        if (u > *(u_neigh + u_deg-1)) {
          high = u_neigh + u_deg;
        }
        else {
          high = u_neigh + u_deg-1;
        }
      }
      while (high>low) {
        int* pivot = (low + (high-low)/2);
        if (*pivot<u) {
          low = pivot + 1;
        }
        else {
          high = pivot;
        }
      }
      vcandidates_start = high;
    }
    
    // Calculate variable m and k
    double r1 = 0.05;
    double r2 = 0.25;
    double r3 = 1.0;
    int64_t mr;
    int64_t kr = 1;
    if (u_deg < (r1 * max_degree)) {
      mr = (max_degree * tres) * r1;
    }
    else if (u_deg < (r2 * max_degree)) {
      mr = (max_degree * tres) * r2;
    }
    else {
      mr = (max_degree * tres) * r3;
    }
    if (max_degree!=0) {
      kr = (mr / max_degree) * 0.69314718056;
    }
    
    BloomSet<NodeID> * neigh = new BloomSet<NodeID>(g.out_neigh(u).begin(), vcandidates_start, mr, kr, seeds, g.num_nodes());
    
    sets.push_back(neigh);
    hstart.push_back(vcandidates_start);
    // cout << "hstart size: " << hstart.size() << endl;
  }
  return hstart;
 
}

template <typename NodeID>
void addSetsStdSet(vector<StdSet<NodeID>*> &sets, const Graph &g)
{
  for (NodeID u = 0; u < g.num_nodes(); ++u) {
    StdSet<NodeID> *neigh = new StdSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end());
    sets.push_back(neigh);
  }
}

size_t OrderedCount(const Graph &g, float tres, int k_, int m_) {
  Timer t, intersection_t;
  t.Start();
  vector<BloomSet<NodeID>*> setsBF;
  vector<StdSet<NodeID>*> sets;
  vector<NodeID*> hstart;
  hstart = addSets(setsBF, g, tres, k_, m_);
  addSetsStdSet(sets, g);
  t.Stop();
  PrintTime("Datastructure building", t.Seconds());
  
  size_t size = 0;
  for (int i = 0; i < g.num_nodes(); i++) {
    size += sets[i]->size();
  }
  PrintTime("Datastructure size", size);

  size_t total = 0;
  size_t n;
  t.Start();
  #pragma omp parallel reduction(+:total)
  {
  #pragma omp for schedule(dynamic)
    for (NodeID u = 0; u < g.num_nodes(); ++u) {
      for (NodeID* vp = g.out_neigh(u).begin(); vp != g.out_neigh(u).end(); vp++) {
        if (*vp > u) break;
        if ((*setsBF[u]).m != (*setsBF[*vp]).m) {
          n = (*sets[u]).intersect_count(*sets[*vp]);
          //cout << "stdset u: " << u << " v: " << *vp << " intersect_count: " << n << endl;
        }
        else {
          n = (*setsBF[u]).intersect_count(*setsBF[*vp]);
          //cout << "BFset u: " << u << " v: " << *vp << " intersect_count: " << n << endl;
        }
        total += n;
      }
    }
    /*
    for (NodeID u = 0; u < g.num_nodes(); ++u) {
      for (NodeID* vp = g.out_neigh(u).begin(); vp != g.out_neigh(u).end() ; vp++) {
      //int64_t deg = g.out_degree(u);
      //for () {
        int n = (*sets[u]).intersect_count(*sets[*vp]);
        cout << "u: " << u << " intersect_count: " << n << endl;
        total += n;
      }
    }
    */
  }
  total /= 3; // TODO: Fix this bug in other TC variations
  t.Stop();
  PrintTime("Intersection time", t.Seconds());
    for(auto set:sets) {
    delete set;
  }
  std::cout << "triangles: " << total << std::endl;
  return total;
}


// heuristic to see if sufficently dense power-law graph
// bool WorthRelabelling(const Graph &g) {
//   int64_t average_degree = g.num_edges() / g.num_nodes();
//   if (average_degree < 10)
//     return false;
//   SourcePicker<Graph> sp(g);
//   int64_t num_samples = min(int64_t(1000), g.num_nodes());
//   int64_t sample_total = 0;
//   pvector<int64_t> samples(num_samples);
//   for (int64_t trial=0; trial < num_samples; trial++) {
//     samples[trial] = g.out_degree(sp.PickNext());
//     sample_total += samples[trial];
//   }
//   sort(samples.begin(), samples.end());
//   double sample_average = static_cast<double>(sample_total) / num_samples;
//   double sample_median = samples[num_samples/2];
//   return sample_average / 1.3 > sample_median;
// }


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
bool TCVerifier(const Graph &g, size_t approx_total, const CLApp &cli) {
  size_t exact_total = 0;
  #pragma omp parallel for reduction(+ : exact_total) schedule(dynamic, 64)
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
          exact_total++;
      }
    }
  }
  cout << "accuracy: " << (float)((float)approx_total-exact_total)/exact_total*100 << "\%" << endl;
  cout << "approx total: " << approx_total << endl;
  cout << "exact total: " << exact_total << endl;
  if (exact_total != approx_total)
    cout << exact_total << " != " << approx_total << endl;
  return exact_total == approx_total;
}



int main(int argc, char* argv[]) {
  CLSIMDApp cli(argc, argv, "triangle count");
  if (!cli.ParseArgs())
    return -1;
  // cli.set_verify();
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  auto ct = [&cli] (const Graph& g) {
    return OrderedCount(g, cli.treshold(), cli.k(), cli.m());
  };
  BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);
  return 0;
}

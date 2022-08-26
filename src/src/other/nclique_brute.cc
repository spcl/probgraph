// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>

#include "../benchmark.h"
#include "../builder.h"
#include "../command_line.h"
#include "../graph.h"
#include "../pvector.h"


/*
Counts the number of 4-cliques
*/


#ifndef COUNT
  #define COUNT 0
#endif

using namespace std;

// Adapted from https://en.wikibooks.org/wiki/Algorithm_Implementation/Search/Binary_search#C++
// Uses a logarithmic number of comparisons, but std::distance still needs
// a linear number of increments
bool is_connected(NodeID u, NodeID v, const Graph &g, size_t *counter=nullptr) {
  NodeID lo;
  NodeID hi;
  // Iterate over the node with lower degree
  #if COUNT
    *counter = *counter + 1;
  #endif
  if (g.out_degree(u) > g.out_degree(v)) {
    lo = v; hi = u;
  } else {
    lo = u; hi = v;
  }

  auto begin = g.out_neigh(lo).begin();
  auto end = g.out_neigh(lo).end();

  // Keep halving the search space until we reach the end of the vector
  while(begin < end) {
    #if COUNT
      *counter = *counter + 1;
    #endif
    // Find the median value between the iterators
    auto middle = begin + (std::distance(begin, end) / 2);

    // Re-adjust the iterators based on the median value
    #if COUNT
      *counter = *counter + 1;
    #endif
    if (*middle == hi) {
      return true;
    }
    else if (*middle > hi) {
      end = middle;
    }
    else {
      begin = middle + 1;
    }
    #if COUNT // for the else if comparison
      *counter = *counter + 1;
    #endif
  }
  return false;
}

// Count a 4-clique only if v1 < v2 < v3 < v4
size_t OrderedCount(const Graph &g) {
  // Instrumentation------------------------------------//
  size_t count1 = 0; // main comparison counter         //
  size_t count2 = 0; // is_connected comparison counter //
  Timer t1; // total timer                              //
  Timer t2; // is_connected timer                       //
  float t2_time = 0;                                    //
  //----------------------------------------------------//

  // size_t count1 = 0; // main comparison counter
  // size_t count2 = 0; // is_connected comparison counter

  size_t total = 0;
  t1.Start();
  #if !COUNT
    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
  #endif
  for (NodeID v1=0; v1 < g.num_nodes(); v1++) {
    for (auto v2 : g.out_neigh(v1)) {
      #if COUNT
        count1++;
      #endif
      if (v2 > v1) {
        for (auto v3 : g.out_neigh(v2)) {
          #if COUNT
            count1++;
          #endif
          if (v3 > v2) {
            for (auto v4 : g.out_neigh(v3)) {
              #if COUNT
                count1++;
              #endif
              if (v4 > v3) {
                // Only three out of six edges need to be checked at this point
                t2.Start();
                bool is_clique = is_connected(v1, v3, g, &count2) &&
                                 is_connected(v1, v4, g, &count2) &&
                                 is_connected(v2, v4, g, &count2);
                t2.Stop();
                t2_time += t2.Seconds();
                if (is_clique) {
                  total += 1;
                }
              }
            }
          }
        }
      }
    }
  }
  t1.Stop();
  cout << "total 4-cliques: " << total << endl;
  PrintTime("total loop time", t1.Seconds());
  cout << "total is_connected time: " << t2_time << endl;
  #if COUNT
    cout << "total comparisons: " << count1 + count2 << endl;
    cout << "main comparisons: " << count1 << endl;
    cout << "is_connected comparisons: " << count2 << endl;
  #endif
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
size_t Hybrid(const Graph &g) {
  if (WorthRelabelling(g))
    return OrderedCount(Builder::RelabelByDegree(g));
  else
    return OrderedCount(g);
}


void PrintTriangleStats(const Graph &g, size_t total_4cliques) {
  cout << total_4cliques << " 4-cliques" << endl;
}


bool Verifier(const Graph &g, size_t test_total, const CLApp &cli) {
  return true;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "4-clique count");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  BenchmarkKernel(cli, g, Hybrid, PrintTriangleStats, Verifier);
  return 0;
}

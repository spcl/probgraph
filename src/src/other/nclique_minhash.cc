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
Counts the number of 4-cliques
*/

#ifndef COUNT
  #define COUNT 0
#endif

using namespace std;

void addSets(vector<KMinHashSet<NodeID>*> &sets, const Graph &g, float threshold){
  struct timeval time;
  gettimeofday(&time, NULL);
  int seed = time.tv_usec;
  int64_t max_degree = 0;
  // Find max degree and use as value for k in combination with threshold
  for(NodeID u=0; u<g.num_nodes(); ++u) {
    max_degree = std::max(g.out_degree(u), max_degree);
  }
  cout << "max deg: " << max_degree << endl;
  float k = max_degree * threshold;
  cout << "k: " << k << endl;

  for (NodeID u=0; u<g.num_nodes(); u++) {
    std::vector<int> vec (g.out_neigh(u).begin(), g.out_neigh(u).end());
    vec.push_back(u);
    int* begin = &vec[0];
    int* end = &vec[0] + vec.size();
    KMinHashSet<NodeID> * neigh = new KMinHashSet<NodeID>(begin, end, k, seed);
    sets[u] = neigh;
  }
}

// Return the max size of a list of vectors
inline size_t max_vector_size(vector<KMinHashSet<int>*> &sets) {
  size_t max_degree = 0;
  for (auto set_p : sets) {
    size_t size = set_p->size();
    max_degree = max(size, max_degree);
  }
  return max_degree;
}

size_t intersect_count(vector<int> &nodes, vector<KMinHashSet<NodeID>*> &sets, size_t *counter) {
  vector<KMinHashSet<NodeID>*> isect_sets;
  isect_sets.reserve(nodes.size()); 
  for (auto node : nodes) {
    isect_sets.push_back(sets[node]);
  }

  KMinHashSet<NodeID> *set = isect_sets.back();
  isect_sets.pop_back();
  return set->intersect_count_multiple(isect_sets, counter);
}

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

size_t OrderedCount(const Graph &g, float tres) {
  // g.PrintTopology();
  Timer t;
  t.Start();
  vector<KMinHashSet<NodeID>*> sets;
  sets.reserve(g.num_nodes());
  addSets(sets, g, tres);
  t.Stop();
  PrintTime("Datastructure building", t.Seconds());

  size_t size = 0;
  for(int i = 0; i<g.num_nodes(); i++){
      size+= sets[i]->size();
  }
  PrintTime("Datastructure size", size);

  // Instrumentation------------------------------------//
  size_t count1 = 0; // main comparison counter         //
  size_t count2 = 0; // is_connected comparison counter //
  size_t count3 = 0; // intersection comparison counter //
  Timer t1; // total timer                              //
  Timer t2; // is_connected timer                       //
  Timer t3; // intersect_count timer                    //
  float t2_time = 0;                                    //
  float t3_time = 0;                                    //
  //----------------------------------------------------//

  size_t n_clique = 4;
  size_t total = 0;
  t1.Start();
  #if !COUNT
    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
  #endif
  for (NodeID v1=0; v1 < g.num_nodes(); v1++) {
    bool cond = false;
    for (NodeID v2 : g.out_neigh(v1)) {
      #if COUNT
        count1++;
      #endif
      if (v2 > v1) {
        vector<NodeID> sets1;
        sets1.push_back(v1);
        sets1.push_back(v2);
        t3.Start();
        size_t count = intersect_count(sets1, sets, &count3);
        t3.Stop();
        t3_time += t3.Seconds();
        #if COUNT
          count1++;
        #endif
        if (count >= n_clique) {
          for (NodeID v3 : g.out_neigh(v1)) {
            #if COUNT
              count1++;
            #endif
            if (v3 > v2) {
              t2.Start();
              cond = is_connected(v3, v2, g, &count2);
              t2.Stop();
              t2_time += t2.Seconds();
            }
            #if COUNT
              count1++;
            #endif
            if (v3 > v2 && cond) {
              vector<NodeID> sets2;
              sets2.push_back(v1);
              sets2.push_back(v2);
              sets2.push_back(v3);
              t3.Start();
              count = intersect_count(sets2, sets, &count3);
              t3.Stop();
              t3_time += t3.Seconds();
              #if COUNT
                count1++;
              #endif
              if (count >= n_clique) {
                for (NodeID v4 : g.out_neigh(v1)) {
                  #if COUNT
                    count1++;
                  #endif
                  if (v4 > v3) {
                    t2.Start();
                    cond = is_connected(v4, v3, g, &count2) && is_connected(v4, v2, g, &count2);
                    t2.Stop();
                    t2_time += t2.Seconds();
                  }
                  #if COUNT
                    count1++;
                  #endif
                  if (v4 > v3 && cond) {
                    vector<NodeID> sets3;
                    sets3.push_back(v1);
                    sets3.push_back(v2);
                    sets3.push_back(v3);
                    sets3.push_back(v4);
                    t3.Start();
                    count = intersect_count(sets3, sets, &count3);
                    t3.Stop();
                    t3_time += t3.Seconds();
                    #if COUNT
                      count1++;
                    #endif
                    if (count >= n_clique) {
                      // cout << "count: " << count << endl;
                      // cout << n_clique << "-clique: " << v1 << " " << v2 << " " << v3 << " " << v4 << endl;
                      total += 1;
                    }
                  }
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
  cout << "total intersection_count time: " << t3_time << endl;
  #if COUNT
    cout << "total comparisons: " << count1 + count2 + count3 << endl;
    cout << "main comparisons: " << count1 << endl;
    cout << "is_connected comparisons: " << count2 << endl;
    cout << "intersection comparisons: " << count3 << endl;
  #endif

  for(auto set:sets) delete set; 
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
  // cout << total_triangles << " triangles" << endl;
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
  cout << "acc: " << (float)((float)test_total-total)/total*100  << "\% error: true " << total << " counted " << test_total << endl;
  if (total != test_total)
    cout << total << " != " << test_total << endl;
  return total == test_total;
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
  auto ct = [&cli] (const Graph& g){
     return OrderedCount(g, cli.treshold());
   };
  BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);
  return 0;
}
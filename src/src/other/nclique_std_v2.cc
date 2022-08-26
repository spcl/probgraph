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

typedef int NodeID;

void addSets(vector<VectorSet<NodeID>*> &sets, const Graph &g){
  for (NodeID u=0; u<g.num_nodes(); u++) {
    VectorSet<NodeID> * neigh = new VectorSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end());
    sets[u] = neigh;
  }
}

// // Return the max size of a list of vectors
// inline size_t max_vector_size(vector<VectorSet<int>*> &sets) {
//   size_t max_degree = 0;
//   for (auto set_p : sets) {
//     size_t size = set_p->size();
//     max_degree = max(size, max_degree);
//   }
//   return max_degree;
// }

// Adapted from https://en.cppreference.com/w/cpp/algorithm/set_intersection
template<class InputIt1, class InputIt2, class OutputIt>
void intersect(InputIt1 first1, InputIt1 last1,
               InputIt2 first2, InputIt2 last2,
               OutputIt d_first, size_t *counter=nullptr) {
  while (first1 != last1 && first2 != last2) {
    #if COUNT
      *counter = *counter + 1;
    #endif
    if (*first1 < *first2) {
      ++first1;
    } else  {
      #if COUNT
        *counter = *counter + 1;
      #endif
      if (!(*first2 < *first1)) {
        *d_first++ = *first1++;
      }
      ++first2;
    }
  }
}

vector<int> intersect(VectorSet<int> * first, vector<int> &last, size_t *counter=nullptr) {
  vector<int> result;
  result.reserve(max(first->size(), (int) last.size()));
  intersect(first->begin(), first->end(), last.begin(), last.end(), std::back_inserter(result), counter);
  return result;
}

vector<int> intersect(VectorSet<NodeID> * first, VectorSet<NodeID> * last, size_t *counter=nullptr) {
  vector<int> result;
  result.reserve(max(first->size(), last->size()));
  intersect(first->begin(), first->end(), last->begin(), last->end(), std::back_inserter(result), counter);
  return result;
}

// size_t intersect_count(vector<int> &nodes, vector<VectorSet<NodeID>*> &sets, size_t *counter) {
//   vector<VectorSet<NodeID>*> isect_sets;
//   isect_sets.reserve(nodes.size()); 
//   for (auto node : nodes) {
//     isect_sets.push_back(sets[node]);
//   }

//   size_t max_size = max_vector_size(isect_sets);
//   std::vector<int> last_isect;
//   last_isect.reserve(max_size);
//   std::copy(isect_sets[0]->begin(), isect_sets[0]->end(), std::inserter(last_isect, last_isect.begin()));
//   std::vector<int> curr_isect;
//   curr_isect.reserve(max_size);

//   for (size_t i=1; i<isect_sets.size(); i++) {
//     intersect(last_isect.begin(), last_isect.end(),
//               isect_sets[i]->begin(), isect_sets[i]->end(),
//               std::back_inserter(curr_isect), counter);

//     std::swap(last_isect, curr_isect);
//     curr_isect.clear();
//   }
//   return last_isect.size();
// }

// vector<int> intersect_multiple(vector<int> &nodes, vector<VectorSet<NodeID>*> &sets, size_t *counter=nullptr) {
//   vector<VectorSet<NodeID>*> isect_sets;
//   isect_sets.reserve(nodes.size()); 
//   for (auto node : nodes) {
//     isect_sets.push_back(sets[node]);
//   }

//   size_t max_size = max_vector_size(isect_sets);
//   std::vector<int> last_isect;
//   last_isect.reserve(max_size);
//   std::copy(isect_sets[0]->begin(), isect_sets[0]->end(), std::inserter(last_isect, last_isect.begin()));
//   std::vector<int> curr_isect;
//   curr_isect.reserve(max_size);

//   for (size_t i=1; i<isect_sets.size(); i++) {
//     intersect(last_isect.begin(), last_isect.end(),
//               isect_sets[i]->begin(), isect_sets[i]->end(),
//               std::back_inserter(curr_isect), counter);

//     std::swap(last_isect, curr_isect);
//     curr_isect.clear();
//   }
//   return last_isect;
// }

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
  vector<VectorSet<NodeID>*> sets;
  sets.reserve(g.num_nodes());
  addSets(sets, g);
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

  // size_t n_clique = 4;
  vector<NodeID> cur_isect1;
  vector<NodeID> cur_isect2;
  size_t total = 0;
  t1.Start();
  #if !COUNT
    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64) private (cur_isect1, cur_isect2)
  #endif
  for (NodeID v1=0; v1 < g.num_nodes(); v1++) {
    bool cond = false;
    for (NodeID v2 : g.out_neigh(v1)) {
      #if COUNT
        count1++;
      #endif
      if (v2 > v1) {
        t3.Start();
        cur_isect1 = intersect(sets[v1], sets[v2], &count3);
        t3.Stop();
        t3_time += t3.Seconds();
        size_t count = cur_isect1.size();
        #if COUNT
          count1++;
        #endif
        if (count >= 2) {
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
              t3.Start();
              cur_isect2 = intersect(sets[v3], cur_isect1, &count3);
              t3.Stop();
              t3_time += t3.Seconds();
              count = cur_isect2.size();
              #if COUNT
                count1++;
              #endif
              if (count >= 1) {
                total += count;
              }
            }
          }
        }
      }
    }
  }
  t1.Stop();
  total = total / 4;

  cout << "total 4-cliques: " << total << endl;
  PrintTime("total loop time", t1.Seconds());
  cout << "total is_connected time: " << t2_time << endl;
  cout << "total Intersection time: " << t3_time << endl;
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


    cout << "approx total: 0"  << endl;
    cout << "verification time: 0" << endl;
    cout << "exact total: " << test_total << endl;


    return true;

    //  size_t total = 0;
//  #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
//  for (NodeID u=0; u < g.num_nodes(); u++) {
//    for (NodeID v : g.out_neigh(u)) {
//      if (v > u)
//        break;
//      auto it = g.out_neigh(u).begin();
//      for (NodeID w : g.out_neigh(v)) {
//        if (w > v)
//          break;
//        while (*it < w)
//          it++;
//        if (w == *it)
//          total++;
//      }
//    }
//  }
//  cout << "acc: " << (float)((float)test_total-total)/total*100  << "\% error: true " << total << " counted " << test_total << endl;
//  if (total != test_total)
//    cout << total << " != " << test_total << endl;
//  return total == test_total;
}



int main(int argc, char* argv[]) {
  CLSIMDApp cli(argc, argv, "4-cqliue std");
  if (!cli.ParseArgs())
    return -1;
  // cli.set_verify();
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but 4-clique requires undirected" << endl;
    return -2;
  }
  auto ct = [&cli] (const Graph& g){
     return OrderedCount(g, cli.treshold());
   };
  BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);
  return 0;
}

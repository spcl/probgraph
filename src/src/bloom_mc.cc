// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "bloom_filter.hpp"
#include <boost/heap/binomial_heap.hpp>
#include <boost/optional.hpp>

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

typedef vector<NodeID> set;
typedef std::pair<NodeID, NodeID> bin_node;
typedef boost::heap::binomial_heap<bin_node>::handle_type handle_t;
vector<bloom_filter> generate_filters(const double acc, const Graph &g){
  const NodeID n = g.num_nodes();
  int max_deg = -1;
  for(auto u = 0; u<n; u++){
    const auto deg_u = g.out_degree(u);
    if(deg_u > max_deg) max_deg=deg_u;
  }
  vector<bloom_filter> filters;
  filters.reserve(n);
  auto seed = random_device{}();
  for(NodeID u = 0; u<n; u++){
    bloom_parameters param;
    param.projected_element_count=max_deg;
    param.false_positive_probability=acc;
    param.random_seed=seed;
    param.compute_optimal_parameters();
    filters.push_back(bloom_filter(param));
    (filters[u]).insert(g.out_neigh(u).begin(), g.out_neigh(u).end());
  }
  bool cond = filters[(int32_t)30].contains((int32_t)0);
  return filters;
}

set degeneracyOrdering(const Graph &g){
    boost::heap::binomial_heap<bin_node> L;
    vector<boost::optional<handle_t>> handles(g.num_nodes());
    set res = set();

    for (NodeID n = 0; n<g.num_nodes(); n++){
        bin_node no = bin_node(-g.out_degree(n),n);
        handles[n] = L.push(no);
    }
    while(!L.empty()){
        bin_node u = L.top();
        L.pop();
        handles[u.second] = boost::none;
        res.push_back(u.second);
        for(auto v : g.out_neigh(u.second))
            if(handles[v]){
            (*(handles[v].get())).first++;
            L.increase(handles[v].get());
            }
    }
    reverse(res.begin(), res.end()); 
    return res;
}

set BronKerboschPivot(set& R,set& P, set& X, const vector<bloom_filter> &filters, const Graph &g) {
    if (P.empty() && X.empty()) return R;
    NodeID u;
    if (!P.empty()){
        u = P.back();
    }else{
        u = X.back();
    }
    set ret;
    for(auto v = P.rbegin(); v!= P.rend(); v++){
        bool cond = filters[u].contains(v);
        if (cond){
            continue;
        } 
        //intersect P,N((*v))
        set Pbar = set();
        for(auto k = P.begin(); k!=P.end(); k++){
            cond = filters[(*v)].contains(*k);
            if(cond)
                Pbar.push_back(*k);
        }
        set Rbar = set(R);
        Rbar.push_back((*v));
        //intersect (X,N((*v)))
        set Xbar = set();
        for(auto k = X.begin(); k!=X.end(); k++){
            if(filters[(*v)].contains(*k))
                Xbar.push_back(*k);
        }       
        ret = BronKerboschPivot(Rbar, Pbar, Xbar, filters, g);
        
        if (ret.size()>0){
            break;
        }
        X.push_back((*v));
        P.pop_back();
    }
    return ret;
}
 
set intersection(NodeID vi,set::iterator  begin, set::iterator end, const vector<bloom_filter> filters){
    set ret;
    for(auto i = begin;i!=end; i++){
        if(filters[vi].contains(*i)){
            ret.push_back(*i);
        }
    }
    return ret;
}

vector<NodeID> MaximumClique(const Graph &g, const vector<bloom_filter> &filters) {
    set ordering = degeneracyOrdering(g);
    set::iterator vi = ordering.begin();
    // for v0 P is N(v0), X is empty, R is {v0}
    set P = set(g.out_neigh(*vi).begin(), g.out_neigh(*vi).end());
    set X = set();
    set R = set();
    R.push_back(*vi);
    return BronKerboschPivot(R,P,X, filters, g);
    vi = vi + 1;

    // middle part
    for(;vi!=ordering.end()-1; vi++){
      set P = set();
      set X = set();
      set R = set();
      R.push_back(*vi);

      P = intersection(*vi, vi+1, ordering.end(), filters);
      X = intersection(*vi, ordering.begin(), vi-1, filters);
      return BronKerboschPivot(R,P,X,filters, g);
    }     

    // last for v(n-1), P is empty, X is N(v(n-1)), R={v(n-1)}
    set Plast = set();
    set Xlast = set(g.out_neigh(*vi).begin(),g.out_neigh(*vi).end());
    set Rlast = set();
    Rlast.push_back(*vi);
    return BronKerboschPivot(Rlast,Plast,Xlast,filters, g); 
}

void PrintMaxCliqueStats(const Graph &g, set largest_clique) {
  cout << "largest clique size: "<< largest_clique.size() << " elements:";
  for(NodeID i:largest_clique){
    cout<<" " <<i;
  }
  cout << endl;
}

// always returns true at the moment
bool MCVerifier(const Graph &g, set mc) {
  return true;
//  size_t total = 0;
//  vector<NodeID> intersection;
//  intersection.reserve(g.num_nodes());
//  for (NodeID u : g.vertices()) {
//    for (NodeID v : g.out_neigh(u)) {
//      auto new_end = set_intersection(g.out_neigh(u).begin(),
//                                      g.out_neigh(u).end(),
//                                      g.out_neigh(v).begin(),
//                                      g.out_neigh(v).end(),
//                                      intersection.begin());
//      intersection.resize(new_end - intersection.begin());
//      total += intersection.size();
//    }
//  }
//  total = total / 6;  // each triangle was counted 6 times
//  if (total != test_total)
//    cout << total << " != " << test_total << " " << (float)((float)test_total-total)/total*100 << "\% error" << endl;
//  return total == test_total;
}

int main(int argc, char* argv[]) {	
  CLBloom_tc cli(argc, argv, "bloom triangle count", 1e-3, 1);
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  cout << "Accuracy: " << (float) cli.accuracy() << endl;
  const auto filters = generate_filters(cli.accuracy(), g);
  auto mc = [&filters] (const Graph &g) {
    return MaximumClique(g, filters);
  };
  BenchmarkKernel(cli, g, mc, PrintMaxCliqueStats, MCVerifier);
  return 0;
}

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
#include "sets.hpp"

/*
Return the number of MAXIMAL CLIQUES.

IMPORTANT: this code is not correct, for example for the input `./mc -g 2 -a`
it returns 1 although the correct number of maximal cliques is 2.

Debug statements can be removed once the code is correct.
*/


using namespace std;

typedef std::pair<NodeID, NodeID> bin_node;
typedef boost::heap::binomial_heap<bin_node>::handle_type handle_t;


vector<NodeID> DegeneracyOrdering(const Graph &g) {
  boost::heap::binomial_heap<bin_node> L;
  vector<boost::optional<handle_t>> handles(g.num_nodes());
  vector<NodeID> res;
  res.reserve(g.num_nodes());

  for (NodeID n = 0; n<g.num_nodes(); n++) {
    bin_node no = bin_node(-g.out_degree(n),n);
    handles[n] = L.push(no);
  }
  while (!L.empty()) {
    bin_node u = L.top();
    L.pop();
    handles[u.second] = boost::none;
    res.push_back(u.second);
    for(auto v : g.out_neigh(u.second)) {
      if (handles[v]) {
        (*(handles[v].get())).first++;
        L.increase(handles[v].get());
      }
    }
  }
  reverse(res.begin(), res.end()); 
  return res;
}

StdSet<NodeID> GetNeighbors(const NodeID u, const Graph &g) {
  StdSet<NodeID> result;
  for (NodeID neighbor : g.out_neigh(u)) {
    result.insert(neighbor);
  }
  return result;
}

void PrintArgs(StdSet<NodeID>& R, StdSet<NodeID>& P, StdSet<NodeID>& X) {
  cout << "R: ";
    for (auto it = R.begin(); it != R.end(); it++) {
      cout << *it << " ";
    }
    cout << ";  ";
    cout << "P: ";
    for (auto it = P.begin(); it != P.end(); it++) {
      cout << *it << " ";
    }
    cout << ";  ";
    cout << "X: ";
    for (auto it = X.begin(); it != X.end(); it++) {
      cout << *it << " ";
    }
    cout << endl;
}

size_t BronKerboschPivot(StdSet<NodeID>& R, StdSet<NodeID>& P, StdSet<NodeID>& X, const Graph &g) {
  cout << "-----> entering BronKerboschPivot with " << endl;
  PrintArgs(R, P, X);
  if (P.empty() && X.empty()) {
    cout << "!FOUND A MAXIMAL CLIQUE!" << endl;
    return 1;
  }
  // Choose a pivot
  NodeID pivot;
  if (!P.empty()) {
    pivot = *P.begin();
  } else {
    pivot = *X.begin();
  }
  cout << "pivot is " << pivot << endl;

  size_t num_cliques = 0;
  // Iterate through P \ N(u)
  StdSet<NodeID> pivot_neighbors = GetNeighbors(pivot, g);

  if (P.empty()) { cout << "P is empty, skipping BKP loop" << endl; }
  for (auto it = P.begin(); it != P.end(); it++) {
    NodeID v = *it;
    cout << "loop start with v = " << v << endl;
    if (!pivot_neighbors.contains(v)) {
      StdSet<NodeID> v_neighbors = GetNeighbors(v, g);
      StdSet<NodeID> P_prime = P.intersect(v_neighbors);
      StdSet<NodeID> R_prime = StdSet<NodeID>(R.begin(), R.end());
      R_prime.insert(v);
      StdSet<NodeID> X_prime = X.intersect(v_neighbors);

      num_cliques += BronKerboschPivot(R_prime, P_prime, X_prime, g);

      // Remove v from P
      it = P.erase(it);
      X.insert(v);
      // Prevent infinite loop
      if (it == P.end()) {
        break;
      }
    }
    cout << "loop end" << endl;
  }
  cout << "returning num_cliques" << endl;
  return num_cliques;
}

size_t MaximumClique(const Graph &g) {
  g.PrintTopology();
  StdSet<NodeID> R;
  StdSet<NodeID> P;
  StdSet<NodeID> X;

  // Initially, P incluldes all nodes in g
  for (NodeID u : g.vertices()) {
    P.insert(u);
  }

  size_t total_cliques = 0;

  // For each vertex in a degeneracy ordering of g
  vector<NodeID> ordering = DegeneracyOrdering(g);
  
  cout << "degeneracy ordering is: ";
  for (auto it = ordering.begin(); it != ordering.end(); it++) {
    cout << *it << " ";
  }
  cout << endl;

  for (auto it = ordering.begin(); it != ordering.end(); it++) {
    NodeID u = *it;
    StdSet<NodeID> u_neighbors = GetNeighbors(u, g);
    P = StdSet<NodeID>(it+1, ordering.end());
    P.intersect(u_neighbors);
    for (auto it = P.begin(); it != P.end(); it++) {

    }
    X = StdSet<NodeID>(ordering.begin(), it);
    X.intersect(u_neighbors);
    R.insert(u);

    cout << "=====> calling from MaximumClique with" << endl;
    PrintArgs(R, P, X);
    
    total_cliques += BronKerboschPivot(R, P, X, g);
    R.erase(u);
  }
  cout << "returning total_cliques" << endl;
  return total_cliques;
}

void PrintMaxCliqueStats(const Graph &g, size_t total_cliques) {
  cout << total_cliques << " maximal cliques" << endl;
}

bool MCVerifier(const Graph &g, size_t total_cliques) {
  // TODO
  return true;
}

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "maximal cliques");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but mc requires undirected" << endl;
    return -2;
  }
  BenchmarkKernel(cli, g, MaximumClique, PrintMaxCliqueStats, MCVerifier);
  return 0;
}

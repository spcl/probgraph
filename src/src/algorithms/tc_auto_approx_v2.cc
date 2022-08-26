#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "../benchmark.h"
#include "../builder.h"
#include "../command_line.h"
#include "../graph.h"
#include "../pvector.h"

#include "../sisa.h"
#include <omp.h>

/*
A modified version of the following source code:

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

size_t OrderedCount(const Graph &g, std::string graphName, int threads, float precision) {
  double tc_time = -1;
  Timer t;
  t.Start();
  size_t total = 0;
  int processing_threshold = (int)((float)g.num_nodes() * precision);


  int nthreads = omp_get_max_threads();
  //std::vector<NodeID> messages[g.num_nodes()][nthreads];
  std::vector<NodeID> messages[g.num_nodes()];
  std::cout << "nthreads: " << nthreads << std::endl;
  // super step 0
#pragma omp parallel for schedule(dynamic, 64)
  for (NodeID u = 0; u < g.num_nodes(); u++) {
    for (NodeID v : g.out_neigh(u)) {
      //if (u < v) {
      if (u > v) {
        messages[u].push_back(v);
      }
    }
    messages[u].shrink_to_fit();
  }

  //std::vector<NodeID> messages_s1[g.num_nodes()][nthreads];
  std::unordered_map<NodeID, std::vector<NodeID> > messages_s1[nthreads];
  // super step 1
#pragma omp parallel for schedule(dynamic, 64)
  for (NodeID u = 0; u < g.num_nodes(); u++) {
    int tid = omp_get_thread_num();
    for (NodeID m = 0; m < messages[u].size(); m += precision) {
      for (NodeID v : g.out_neigh(u)) {
        if (messages[u][m] < u && u < v) {
          //messages_s1[v][tid].push_back(messages[u][m]);
          messages_s1[tid][v].push_back(messages[u][m]);
        }
      }
    }
  }

  // super step 2
#pragma omp parallel for reduction(+:total) schedule(dynamic, 64)
  for (NodeID u = 0; u < g.num_nodes(); u++) {
    NodeID m = 0;
    for (int t = 0; t < nthreads; t++) {
      if (messages_s1[t].find(u) != messages_s1[t].end()) {
        for (; m < messages_s1[t][u].size(); m += precision) {
          for (NodeID v : g.out_neigh(u)) {
            if (messages_s1[t][u][m] == v) {
              total++;
            }
          }
        }
        m = m - messages_s1[t][u].size();
      }
    }
  }
 
  total = total * precision;
  t.Stop();

  PrintTime("BSP TC time", t.Seconds());
  tc_time = t.Seconds();

  auto threshold = precision;
  auto k = 0;
  auto m = 0;
  auto pp_time = 0;
  auto approx_str_size = 0;
  auto initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

  // RRR - this means that a given line is dedicated to the runtime results
  // the columns are as follows:
  // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

  std::cout << "RRR TC BSP2 TC_BSP2 " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " "  << m << " " << pp_time << " " << tc_time << " " << pp_time + tc_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

  std::cout << "SSS TC BSP2 TC_BSP2 " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;


  size_t test_total = total;
  total = 0;
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
  cout << "acc: " << (float)((float)test_total-total)/total*100  << "\% error: true " << total << " counted " << test_total << endl;

  return total;
}


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
    cout << total_triangles << " triangles" << endl;
}

bool TCVerifier_null(const Graph &g, size_t approx_total, const CLApp &cli) {
    return true;
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
    cout << "acc: " << (float)((float)test_total-total)/total*100  << "\% error: true " << total << " counted " << test_total << endl;
    return total == test_total;
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

    auto ct = [&cli] (const Graph& g){
        return OrderedCount(g, cli.getGraphBasename(), cli.getThreadNum(), cli.treshold());
    };

    BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier_null);
    return 0;
}

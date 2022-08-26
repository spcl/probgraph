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

size_t OrderedCount(const Graph &g, float tau, std::string graphName, int threads) {
    double tc_time = -1;
    Timer t;
    t.Start();
    size_t total = 0;
    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
    for (NodeID u=0; u < g.num_nodes(); u++) {
        for (NodeID v : g.out_neigh(u)) {
            if (v > u)
                break;

            auto it1 = g.out_neigh(v).begin();
            auto it2 = g.out_neigh(u).begin();
            auto end1 = g.out_neigh(v).end();
            auto end2 = g.out_neigh(u).end();
            int64_t total_buf = 0;

            while (it1 != end1 && it2 != end2) {
                if (*it1 < *it2) {
                    it1++;
                } else if (*it1 > *it2) {
                    it2++;
                } else {
                    it1++;
                    it2++;
                    total_buf++;
                }
            }

            if (total_buf * 1.0 / std::min(g.out_degree(u), g.out_degree(v)) > tau) {
                total++;
            }
        }
    }

    t.Stop();
    //  PrintTime("TC runtime", t.Seconds());
    PrintTime("Intersection time", t.Seconds());
    tc_time = t.Seconds();

    auto threshold = 0;
    auto k = 0;
    auto m = 0;
    auto pp_time = 0;
    auto approx_str_size = 0;
    auto initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR JP-OV BASE JP-OV_BASE " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " "  << m << " " << tau << " " << pp_time << " " << tc_time << " " << pp_time + tc_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS JP-OV BASE JP-OV_BASE " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << tau << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;


    return total;
}


void PrintClusterStats(const Graph &g, int64_t size) {
    cout << size << " edges in the clustering " <<endl;
}


// Compares with simple serial implementation that uses std::set_intersection
bool TCVerifier(const Graph &g, size_t test_total, const CLApp &cli) {
    return true;
}


int main(int argc, char* argv[]) {
    CLApp cli(argc, argv, "triangle count");
    if (!cli.ParseArgs())
        return -1;
    Builder b(cli);
    Graph g = b.MakeGraph();
    if (g.directed()) {
        cout << "Input graph is directed but tc requires undirected" << endl;
        return -2;
    }

    auto ct = [&cli] (const Graph& g){
        return OrderedCount(g, cli.tau(), cli.getGraphBasename(), cli.getThreadNum());
    };

    BenchmarkKernel(cli, g, ct, PrintClusterStats, TCVerifier);
    return 0;
}

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
#include <unordered_set>


/*
Counts the number of 4-cliques
*/


using namespace std;

// Adapted from https://en.wikibooks.org/wiki/Algorithm_Implementation/Search/Binary_search#C++
// Uses a logarithmic number of comparisons, but std::distance still needs
// a linear number of increments
bool is_connected(NodeID u, NodeID v, const Graph &g, size_t *counter=nullptr) {
    NodeID lo;
    NodeID hi;
    // Iterate over the node with lower degree
    if (g.out_degree(u) > g.out_degree(v)) {
        lo = v; hi = u;
    } else {
        lo = u; hi = v;
    }

    auto begin = g.out_neigh(lo).begin();
    auto end = g.out_neigh(lo).end();

    // Keep halving the search space until we reach the end of the vector
    while(begin < end) {
        // Find the median value between the iterators
        auto middle = begin + (std::distance(begin, end) / 2);

        // Re-adjust the iterators based on the median value
        if (*middle == hi) {
            return true;
        }
        else if (*middle > hi) {
            end = middle;
        }
        else {
            begin = middle + 1;
        }
    }
    return false;
}

// Count a 4-clique only if v1 < v2 < v3 < v4
size_t OrderedCount(const Graph &g, std::string graphName, int threads, float precision) {
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
    std::unordered_set<NodeID> snodes[g.num_nodes()];
    t1.Start();
    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
    for (NodeID v1=0; v1 < g.num_nodes(); v1++) {
        for (NodeID v : g.out_neigh(v1)) {
          snodes[v1].insert(v);
        }
    }

    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
    for (NodeID v1=0; v1 < g.num_nodes(); v1++) {
		for (auto v2 : g.out_neigh(v1)) {
			float rnd_num = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			if(rnd_num <= precision){
				if (v2 > v1) {
					for (auto v3 : g.out_neigh(v2)) {
						if (v3 > v2) {
							for (auto v4 : g.out_neigh(v3)) {
								if (v4 > v3) {
									// Only three out of six edges need to be checked at this point
									/*
									bool is_clique = is_connected(v1, v3, g, &count2) &&
													 is_connected(v1, v4, g, &count2) &&
													 is_connected(v2, v4, g, &count2);
									*/
									bool is_clique = !(snodes[v1].find(v3) == snodes[v1].end());
									if (is_clique) {
									  is_clique = !(snodes[v1].find(v4) == snodes[v1].end());
									  if (is_clique) {
										is_clique = !(snodes[v2].find(v4) == snodes[v2].end());
									  }
									}
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
    }
    t1.Stop();
    cout << "total 4-cliques: " << total << endl;
    PrintTime("total loop time", t1.Seconds());
    cout << "total is_connected time: " << t2_time << endl;
    auto alg_time = t1.Seconds();

    auto threshold = precision;
    auto k = 0;
    auto m = 0;
    auto pp_time = 0;
    auto approx_str_size = 0;
    auto initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR 4C PGP 4C_PGP " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " "  << m << " " << pp_time << " " << alg_time << " " << pp_time + alg_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS 4C PGP 4C_PGP " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;


    return total;
}


void PrintTriangleStats(const Graph &g, size_t total_4cliques) {
    cout << total_4cliques << " 4-cliques" << endl;
}


bool Verifier(const Graph &g, size_t test_total, const CLApp &cli) {
    return true;
}


int main(int argc, char* argv[]) {
    CLSIMDApp cli(argc, argv, "4-clique count");
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
    BenchmarkKernel(cli, g, ct, PrintTriangleStats, Verifier);
    return 0;
}

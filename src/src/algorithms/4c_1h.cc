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
#include <unordered_set>


/*
Counts the number of 4-cliques
*/

using namespace std;

void addSets(vector<OneHashSet<NodeID>*> &sets, const Graph &g, float threshold, int& k_){
    struct timeval time;
    gettimeofday(&time, NULL);
    int seed = time.tv_usec;
    int64_t max_degree = 0;
    // Find max degree and use as value for k in combination with threshold
    for(NodeID u=0; u<g.num_nodes(); ++u) {
        max_degree = std::max(g.out_degree(u), max_degree);
    }

    int k = k_;
    if (k <= 0) k = max_degree * threshold;
    if (k < 1) k = 1;
    k_ = k;
    std::cout << "max deg: " << max_degree << std::endl;
    std::cout << "k: " << k << std::endl;

    for (NodeID u=0; u<g.num_nodes(); u++) {
        std::vector<long int> vec (g.out_neigh(u).begin(), g.out_neigh(u).end());
        vec.push_back(u);
        long int* begin = &vec[0];
        long int* end = &vec[0] + vec.size();
        OneHashSet<NodeID> * neigh = new OneHashSet<NodeID>(begin, end, k, seed);
        sets[u] = neigh;
    }
}


size_t intersect_count(vector<long int> &nodes, vector<OneHashSet<NodeID>*> &sets, size_t *counter) {
    vector<OneHashSet<NodeID>*> isect_sets;
    isect_sets.reserve(nodes.size());
    for (auto node : nodes) {
        isect_sets.push_back(sets[node]);
    }


    OneHashSet<NodeID> *set = isect_sets.back();
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

size_t OrderedCount(const Graph &g, float threshold_, int k_, std::string graphName, int threads) {
    // g.PrintTopology();
    Timer t;
    t.Start();
    vector<OneHashSet<NodeID>*> sets;
    sets.reserve(g.num_nodes());
    float threshold = threshold_;
    int k = k_;
    addSets(sets, g, threshold, k);
    t.Stop();
    PrintTime("Datastructure building", t.Seconds());
    auto pp_time = t.Seconds();

    size_t size = 0;
    for(int i = 0; i<g.num_nodes(); i++){
        size+= sets[i]->size();
    }
    PrintTime("Datastructure size", size);
    double approx_str_size = size / (1024.0 * 1024.0); // MB
    double initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

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
        bool cond = false;
        for (NodeID v2 : g.out_neigh(v1)) {
            if (v2 > v1) {
                vector<NodeID> sets1 = {v1, v2};
                size_t count = intersect_count(sets1, sets, &count3);
                if (count >= n_clique) {
                    for (NodeID v3 : g.out_neigh(v1)) {
                        if (v3 > v2) {
                            cond = !(snodes[v3].find(v2) == snodes[v3].end());
                        }
                        if (v3 > v2 && cond) {
                            vector<NodeID> sets2 = {v1, v2, v3};
                            count = intersect_count(sets2, sets, &count3);
                            if (count >= n_clique) {
                                for (NodeID v4 : g.out_neigh(v1)) {
                                    if (v4 > v3) {
                                        cond = !(snodes[v4].find(v3) == snodes[v4].end());
                                        if (cond) {
                                            cond = !(snodes[v4].find(v2) == snodes[v4].end());
                                        }
                                    }
                                    if (v4 > v3 && cond) {
                                        vector<NodeID> sets3 = {v1, v2, v3, v4};
                                        count = intersect_count(sets3, sets, &count3);
                                        if (count >= n_clique) {
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
    auto alg_time = t1.Seconds();

    int m = 0;

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (1-Hash parameter)] [k (another 1-Hash parameter-number of min hashes] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR 4C 1H 4C_1H " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << pp_time << " " << alg_time << " " << pp_time + alg_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (1-Hash parameter)] [k (another 1-Hash parameter-number of min hashes] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS 4C 1H 4C_1H " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;

    printf("ooo -------------------------------\n");

    for(auto set:sets) delete set;
    return total;
}


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
    CLSIMDApp cli(argc, argv, "4 cliques");
    if (!cli.ParseArgs())
        return -1;

    Builder b(cli);
    Graph g = b.MakeGraph();
    if (g.directed()) {
        cout << "Input graph is directed but 4C requires undirected" << endl;
        return -2;
    }
    auto ct = [&cli] (const Graph& g){
        return OrderedCount(g, cli.treshold(), cli.k(), cli.getGraphBasename(), cli.getThreadNum());
    };
    BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);
    return 0;
}

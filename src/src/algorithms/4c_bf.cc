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

std::vector<NodeID*> addSets(vector<BloomSet<NodeID>*> &sets, const Graph &g, float tres, int& k_, int& m_) {
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
    if (m<32) {
        m = 32;
    }

    if (k <= 0){
        if (max_degree!=0) {
            k = m/max_degree*0.69314718056; // m/n*ln(2)
        }
    }
    if (k<1) {
        k = 1;
    }

    cout << "m: " << m << endl;
    cout << "k: " << k << endl;
    k_ = k;
    m_ = m;
    struct timeval seed;
    gettimeofday(&seed, NULL);
    srand(seed.tv_usec);
    std::vector<uint32_t> seeds;
    for (int i = 0; i<k; i++) {
        seeds.push_back(rand());
    }


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
            long int* low = u_neigh + oldind;
            long int* high;
            if (ind<u_deg) {
                high = u_neigh + ind;
            }
            else {
                if (u > *(u_neigh + u_deg-1)) {
                    high = u_neigh + u_deg;
                } else {
                    high = u_neigh + u_deg-1;
                }
            }
            while (high>low) {
                long int* pivot = (low + (high-low)/2);
                if (*pivot<u) {
                    low = pivot + 1;
                } else {
                    high = pivot;
                }
            }
            vcandidates_start = high;
        }

        BloomSet<NodeID> * neigh = new BloomSet<NodeID>(g.out_neigh(u).begin(), vcandidates_start, m, k, seeds, g.num_nodes());

        sets.push_back(neigh);
        hstart.push_back(vcandidates_start);
    }
    return hstart;
}

size_t intersect_count(vector<long int> &nodes, vector<BloomSet<NodeID>*> &sets, size_t *counter) {
    vector<BloomSet<NodeID>*> isect_sets;
    isect_sets.reserve(nodes.size());
    for (auto node : nodes) {
        isect_sets.push_back(sets[node]);
    }

    BloomSet<NodeID> *set = isect_sets.back();
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

size_t OrderedCount(const Graph &g, float tres, int k_, int m_, std::string graphName, int threads) {
    // g.PrintTopology();
    Timer t;
    t.Start();
    vector<BloomSet<NodeID>*> sets;
    sets.reserve(g.num_nodes());
    int k = k_;
    int m = m_;
    addSets(sets, g, tres, k, m);
    t.Stop();
    PrintTime("Datastructure building", t.Seconds());
    double pp_time = t.Seconds();

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
    Timer t4; // intersect_count timer                    //
    float t2_time = 0;                                    //
    float t3_time = 0;                                    //
    //----------------------------------------------------//

    size_t n_clique = 4;
    size_t total = 0;
    std::unordered_set<NodeID> snodes[g.num_nodes()];
    t1.Start();
    t4.Start();
    #pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
    for (NodeID v1=0; v1 < g.num_nodes(); v1++) {
        for (NodeID v : g.out_neigh(v1)) {
          snodes[v1].insert(v);
        }
    }
    t4.Stop();
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
    PrintTime("unordered set build time", t4.Seconds());
    //cout << "total is_connected time: " << t2_time << endl;
    //cout << "total intersection_count time: " << t3_time << endl;
    double alg_time = t1.Seconds();

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [treshold (BF parameter)] [b (another BF parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR 4C BF 4C_BF " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << tres << " " << k << " " << m << " " << pp_time << " " << alg_time << " " << pp_time + alg_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [treshold (BF parameter)] [b (another BF parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS 4C BF 4C_BF " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << tres << " " << k << " " << m << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;

    printf("ooo -------------------------------\n");

    for(auto set:sets) delete set;
    return total;
}


void PrintCountStats(const Graph &g, size_t total_count) {
    // cout << total_count << "4c" << endl;
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
    // cli.set_verify();
    Builder b(cli);
    Graph g = b.MakeGraph();
    if (g.directed()) {
        cout << "Input graph is directed but 4C requires undirected" << endl;
        return -2;
    }
    auto ct = [&cli] (const Graph& g){
        return OrderedCount(g, cli.treshold(), cli.k(), cli.m(), cli.getGraphBasename(), cli.getThreadNum());
    };
    BenchmarkKernel(cli, g, ct, PrintCountStats, TCVerifier);
    return 0;
}

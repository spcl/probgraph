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
#include <sys/time.h>
#include "../graph.h"
#include "../pvector.h"
#include "../set_operation.hpp"
#include "../sets.hpp"
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
inline NodeID* get_vcandidates_start(const Graph &g, const NodeID u) {
    NodeID* uneigh = g.out_neigh(u).begin();
    NodeID nuneigh = g.out_degree(u);
    NodeID* vcandidates_start;
    if (nuneigh==0) {
        vcandidates_start = uneigh;
    } else {
        int oldind = 0;
        int ind = 0;
        while (ind < nuneigh) {
            if (*(uneigh+ind)>u) {
                break;
            } else {
                oldind = ind;
                ind = ((ind+1)<<1)-1; // 2x(ind+1) - 1
            }
        }
        long int* low = uneigh + oldind;
        long int* high;
        if (ind<nuneigh) {
            high = uneigh + ind;
        }
        else {
            if (u > *(uneigh + nuneigh-1)) {
                high = uneigh + nuneigh;
            } else {
                high = uneigh + nuneigh-1;
            }
        }
        while (high>low) {
            long int* pivot = (low + (high-low)/2);
            if (*pivot<u) low = pivot + 1;
            else high = pivot;
        }
        vcandidates_start = high;
    }
    return vcandidates_start;
}


int64_t OrderedCount(const Graph &g, float treshold, int k_, float tau, std::string graphName, int threads) {
    std::cout << "threshold: " << treshold << std::endl;
    double alg_time = -1;
    double pp_time = -1;
    double approx_str_size = -1;
    double initial_csr_size = -1;
    Timer t;
    t.Start();
    struct timeval time;
    gettimeofday(&time, NULL);
    int seed = time.tv_usec;
    // One KMinHashSet per node (neighborhood)
    vector<KMinHashSet<NodeID>*> sets;
    sets.reserve(g.num_nodes());
    vector<NodeID*> hstart;
    hstart.reserve(g.num_nodes());
    int64_t max_degree = 0;

    // Find max degree and use as value for k in combination with threshold
    for(NodeID u=0; u<g.num_nodes(); ++u) {
        max_degree = std::max(g.out_degree(u), max_degree);
    }

    int k = k_;
    if (k <= 0) k = max_degree * treshold;
    if (k < 1) k = 1;

    std::cout << "max deg: " << max_degree << std::endl;
    std::cout << "k: " << max_degree * treshold << std::endl;

    for(NodeID u = 0; u<g.num_nodes(); ++u) {
        hstart[u] = get_vcandidates_start(g,u);
        KMinHashSet<NodeID> * neigh = new KMinHashSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end(), max_degree * treshold, seed);
        sets[u] = neigh;
    }
    t.Stop();
    PrintTime("Datastructure building", t.Seconds());
    pp_time = t.Seconds();

    size_t size = 0;
    for(int i = 0; i<g.num_nodes(); i++){
        size+= sets[i]->total_size();
    }
    PrintTime("Datastructure size", size);
    approx_str_size = size / (1024.0 * 1024.0); // MB
    initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

    int64_t total = 0;
    
    t.Start();
    #pragma omp parallel for reduction(+:total) schedule(dynamic, 64)
    for (NodeID u = 0; u < g.num_nodes(); ++u) {
        for (NodeID *vp = g.out_neigh(u).begin(); vp < hstart[u]; vp++) {
            int n = (*sets[u]).intersect_count(*sets[*vp]);
            if(n > tau){
                total += 1;
            }

        }
    }
    t.Stop();
    PrintTime("Intersection time", t.Seconds());
    alg_time = t.Seconds();

    auto m = 0;

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR JP-CN KH JP-CN_KH " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << treshold << " " << k << " " << m << " " << tau << " " << pp_time << " " << alg_time << " " << pp_time + alg_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS JP-CN KH JP-CN_KH " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << treshold << " " << k << " " << m << " " << tau << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;


    return total;
}


void PrintClusterStats(const Graph &g, int64_t size) {
    cout << size << " edges in the clustering " <<endl;
}



bool ClusterSizeVerifier(const Graph &g, int64_t approx_clusters, const CLApp &cli) {
    EdgeVector exact_clusters;
#pragma omp parallel for  schedule(dynamic, 64)
    for (NodeID u=0; u < g.num_nodes(); u++) {

        for (NodeID v : g.out_neigh(u)) {
            if (v > u)
                break;

            int comm_neigh = 0;

            auto it = g.out_neigh(u).begin();


            for (NodeID w : g.out_neigh(v)) {
                while (*it < w and it < g.out_neigh(u).end())
                    it++;
                if (w == *it){
                    comm_neigh++;
                    if (comm_neigh>cli.tau())
                        break;
                }
            }

            if (comm_neigh>cli.tau())
                #pragma omp critical
                exact_clusters.push_back({v,u});
        }
    }

    cout << "accuracy: " <<  (float)((float)approx_clusters-exact_clusters.size())/exact_clusters.size()*100 << "\%" << endl;
    cout << "approx cluster size: " << approx_clusters << endl;
    cout << "exact cluster size: " << exact_clusters.size() << endl;
    if (exact_clusters.size()!= approx_clusters)
        cout << exact_clusters.size() << " != " << approx_clusters << endl;
    return exact_clusters.size() == approx_clusters;
}



int main(int argc, char* argv[]) {
    CLSIMDApp cli(argc, argv, "MinHash Clustering");
    if (!cli.ParseArgs())
        return -1;
    Builder b(cli);
    Graph g = b.MakeGraph();
    if (g.directed()) {
        cout << "Input graph is directed but clustering requires undirected" << endl;
        return -2;
    }

    auto clusteringcalc = [&cli] (const Graph& g) {
        return OrderedCount(g, cli.treshold(), cli.k(), cli.tau(), cli.getGraphBasename(), cli.getThreadNum());
    };

    BenchmarkKernel(cli, g, clusteringcalc, PrintClusterStats, ClusterSizeVerifier);
    return 0;
}

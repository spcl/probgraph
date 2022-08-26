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

    k_ = k;
    m_ = m;

    cout << "m: " << m << endl;
    cout << "k: " << k << endl;
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
                }
                else {
                    high = u_neigh + u_deg-1;
                }
            }
            while (high>low) {
                long int* pivot = (low + (high-low)/2);
                if (*pivot<u) {
                    low = pivot + 1;
                }
                else {
                    high = pivot;
                }
            }
            vcandidates_start = high;
        }

        BloomSet<NodeID> * neigh = new BloomSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end(), m, k, seeds, g.num_nodes());

        sets.push_back(neigh);
        hstart.push_back(vcandidates_start);
    }
    return hstart;

}



int64_t OrderedCount(const Graph &g, float tres, int k_, int m_,float tau, std::string graphName, int threads) {
    Timer t;
    t.Start();
    double alg_time = -1;
    double pp_time = -1;
    double approx_str_size = -1;
    double initial_csr_size = -1;
    vector<BloomSet<NodeID>*> sets;
    vector<NodeID*> hstart;
    int k_temp;
    int m_temp;

    k_temp = k_;
    m_temp = m_;
    hstart = addSets(sets, g, tres, k_temp, m_temp);
    t.Stop();
    PrintTime("Datastructure building", t.Seconds());
    pp_time = t.Seconds();

    size_t size = 0;
    for (int i = 0; i<g.num_nodes(); i++) {
        size+= sets[i]->total_size();
    }
    PrintTime("Datastructure size", size);
    approx_str_size = size / (1024.0 * 1024.0); // MB
    initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

    t.Start();
    int64_t total = 0;
    #pragma omp parallel for reduction(+:total) schedule(dynamic, 64)
    for (NodeID u = 0; u < g.num_nodes(); ++u) {
        for (NodeID *vp = g.out_neigh(u).begin(); vp < hstart[u]; vp++) {
            int n = (*sets[u]).intersect_count(*sets[*vp]);
            if(n > tau) {
                total += 1;
            }
        }
    }
    t.Stop();
    PrintTime("Intersection time", t.Seconds());
    alg_time = t.Seconds();

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR JP-CN BF JP-CN-BF " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << tres << " " << k_temp << " " << m_temp << " " << tau << " " << pp_time << " " << alg_time << " " << pp_time + alg_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS JP-CN BF JP-CN-BF " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << tres << " " << k_temp << " " << m_temp << " " << tau << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;

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


    cout << "accuracy: " << 100 -  (float)((float)approx_clusters-exact_clusters.size())/exact_clusters.size()*100 << "\%" << endl;
    cout << "approx cluster size: " << approx_clusters << endl;
    cout << "exact cluster size: " << exact_clusters.size() << endl;
    if (exact_clusters.size()!= approx_clusters)
        cout << exact_clusters.size() << " != " << approx_clusters << endl;
    return exact_clusters.size() == approx_clusters;
}



int main(int argc, char* argv[]) {
    CLSIMDApp cli(argc, argv, "BLOOM Clustering");
    if (!cli.ParseArgs())
        return -1;
//     cli.set_verify();
    Builder b(cli);
    Graph g = b.MakeGraph();
    if (g.directed()) {
        cout << "Input graph is directed but clustering algorithm requires undirected" << endl;
        return -2;
    }
    auto clusteringcalc = [&cli] (const Graph& g) {
        return OrderedCount(g, cli.treshold(), cli.k(), cli.m(), cli.tau(), cli.getGraphBasename(), cli.getThreadNum());
    };

    BenchmarkKernel(cli, g, clusteringcalc, PrintClusterStats, ClusterSizeVerifier);
    return 0;
}

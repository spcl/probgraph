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

#include "../sisa.h"

// pratio and qratio values from minebench
#define pratio 4
#define qratio pratio
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

// TODO: can't use const SISA_Graph as many member functions are non constant!
std::vector<NodeID*> addSets(vector<BloomSet<NodeID>*> &sets, SISA_Graph &g, float tres,  int k_, int m_) {
    //cout << "threshold: " << tres << endl;
    int64_t max_degree = 0;
    for (NodeID u=0; u<g.num_nodes(); ++u) {
        max_degree = std::max(max_degree, g.out_degree(u));
    }
    //cout << "max_degree: " << max_degree << endl;
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

    struct timeval seed;
    gettimeofday(&seed, NULL);
    srand(seed.tv_usec);
    std::vector<uint32_t> seeds;
    for (int i = 0; i<k; i++) {
        seeds.push_back(rand());
    }


    //for(NodeID u = 0; u<g.num_nodes(); u++){
    //    BloomSet<NodeID> * neigh = new BloomSet<NodeID>(g.out_neigh(u).begin(), g.out_neigh(u).end(), m, k);

    //    sets.push_back(neigh);

    //}
    std::vector<NodeID*> hstart; // Stores the end iterator to the range of neighbor vertices
    hstart.reserve(g.num_nodes());
    for (NodeID u = 0; u < g.num_nodes(); ++u) {
        NodeID* u_neigh = g.neigh(u).data();
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
            int* low = u_neigh + oldind;
            int* high;
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
                int* pivot = (low + (high-low)/2);
                if (*pivot<u) {
                    low = pivot + 1;
                }
                else {
                    high = pivot;
                }
            }
            vcandidates_start = high;
        }

        BloomSet<NodeID> * neigh = new BloomSet<NodeID>(g.neigh(u).data(), vcandidates_start, m, k, seeds, g.num_nodes());

        sets.push_back(neigh);
        hstart.push_back(vcandidates_start);
    }
    return hstart;

}



int64_t OrderedCount(const Graph &g, const Graph &cg, float tres, int k_, int m_) {

    SISA_Graph SISAg(g);


    int64_t n = SISAg.num_nodes();
    SSSet E = SISAg.GetEdges();
    int64_t m = E.size();
    int64_t p = m/pratio;
    int64_t q = m/qratio;

    std::random_device rd;
    std::mt19937 eng(0);//rd());
    std::uniform_int_distribution<> distr(1, m/p);
    //how many edges are dropped
    std::vector<bool> visited(m);
    SISA_Graph SISAgprime(cg);

    for (auto e : E) {
        SISAgprime.RemoveEdge(E2U(e,n), E2V(e,n));
    }
    //now SISAgprime G' is a complement graph to G

    auto current = E.begin();
    SSSet discardE;
    for (int64_t i = 0; i < p; i++){ // <- pick p random edges
        int64_t e = distr(eng); // <-random edge offset
        for (int64_t j = 0; j < e; j++){
            ++current;
        }
        discardE.AddVertex(*current);
        SISAgprime.AddEdge(E2U(*current,n), E2V(*current,n));
    }
    //now SISAgprime G' has additionally p random edges from G

    SSSet Er = SISAgprime.GetEdges();
    int64_t mprime = Er.size();
    //std::cout << "\nmprime: " << mprime << ", n: " << n << ", m: " << m << ", p: " << p << "\n";
    assert(mprime == n*(n-1) - m + p);
    std::vector<double> score(q);  //holds similiarity score for q best edges
    std::vector<NodeID> scoredE(q);  //holds q best edges


    vector<BloomSet<NodeID>*> sets;
    vector<NodeID*> hstart;

    hstart = addSets(sets, SISAgprime, tres, k_, m_);


    Timer t,t2;
    float inter_time = 0;
    t.Start();
    t.Stop();
    PrintTime("Datastructure building", t.Seconds());

    size_t size = 0;
    for (int i = 0; i<g.num_nodes(); i++) {
        size+= sets[i]->total_size();
    }
    PrintTime("Datastructure size", size);
    for (auto e : Er) {
        //SimType curScore = S<similarity>(SISAgprime.neigh(E2U(e,n)), SISAgprime.neigh(E2V(e,n)));
        t2.Start();
        int64_t curScore = (*sets[E2U(e, n)]).intersect_count(*sets[E2V(e, n)]);
        t2.Stop();
        inter_time += t2.Seconds();
        int64_t curRank = 0;
        while (curScore > score[curRank]) {
            curRank++;
            if (curRank >= q)
                break;
        }
        for (int i = 0; i < (int)curRank-1; i++){
            score[i] = score[i+1];
            scoredE[i] = scoredE[i+1];
        }
        if (curRank > 0) {
            score[curRank-1] = curScore;
            scoredE[curRank-1] = e;
        }
    }

    std::sort(scoredE.begin(), scoredE.end());
    NodeID *sE = &scoredE[0];
    BloomSet<NodeID> *scoredE_BF = new BloomSet<NodeID>(sE, (sE + scoredE.size()), 256 /* TODO: m? */);
    NodeID *dE = &discardE[0];
    BloomSet<NodeID> *discardE_BF = new BloomSet<NodeID>(dE, (dE + discardE.size()), 256 /* TODO: m? */);
    int64_t eff = (*scoredE_BF).intersect_count(*discardE_BF);
    cout << "eff: " << eff << endl;
    // TODO: fix this
    PrintTime("Intersection time", inter_time);

    return eff;
}


void PrintLinkPredictionStats(const Graph &g, int64_t eff) {
}



bool LinkPredictionVerifier(const Graph &g, int64_t approx_result, const CLApp &cli) {
    cout << "accuracy: 0" << endl;
    cout << "approx total: " << approx_result << endl;
    cout << "exact total: 0" << endl;
    return true;
}


Graph cg;
int main(int argc, char* argv[]) {
    CLSIMDApp cli(argc, argv, "Link Prediction bloom filter");
    if (!cli.ParseArgs())
        return -1;
//     cli.set_verify();
    Builder b(cli);
    Graph g = b.MakeGraph();
    cg = b.MakeClique(g.num_nodes());
    if (g.directed()) {
        cout << "Input graph is directed but link prediction algorithm requires undirected" << endl;
        return -2;
    }
    auto LPcalc = [&cli] (const Graph& g) {
        return OrderedCount(g, cg, cli.treshold(), cli.k(), cli.m());
    };

    BenchmarkKernel(cli, g, LPcalc, PrintLinkPredictionStats, LinkPredictionVerifier);
    return 0;
}

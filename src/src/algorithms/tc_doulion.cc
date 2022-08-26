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

#include "../sisa.h"


using namespace std;

int TCDoulion(const Graph &g, float  p_, std::string graphName, int threads)
{


    SISA_Graph SISAg(g);

    std::random_device rd;
    std::mt19937 eng(22);//rd());

    SSSet E = SISAg.GetEdges();
    int64_t n = SISAg.num_nodes();

    std::uniform_real_distribution<> distr(0.0,1.0);
    float p = p_;

    double pp_time = -1;
    double tc_time = -1;


    Timer t;
    t.Start();

#pragma omp parallel for
    //for(auto& e : E)
    //for (NodeID e = 0; e < SIS; u++) {
    for (size_t i = 0; i < E.size(); i++)
    {
        auto e = E[i];
        if (E2U(e,n) < E2V(e,n))
        {
            if (distr(eng) > p)
            {
                //SISAg.RemoveEdge(E2U(e,n) , E2V(e,n));
                //SISAg.RemoveEdge(E2V(e,n) , E2U(e,n));
                SISAg.markEdgesforRemoval(E2U(e,n), E2V(e,n));
                SISAg.markEdgesforRemoval(E2V(e,n), E2U(e,n));
            }
        }
    }

#pragma omp parallel for schedule(dynamic, 64)
    for (int64_t u = 0; u < n; u++) {
      for (int64_t v = 0; v < SISAg.n_neighs[u]; v++) {
        if (SISAg.td[u][v] == true) {
          SISAg.RemoveEdge(u, SISAg.neigh(u)[v]);
        }
      }
    }

    t.Stop();
    PrintTime("Coloring time", t.Seconds());
    pp_time = t.Seconds();

    //SSSet En = SISAg.GetEdges();
    //std::cout << E.size() << "reduced edge size: " << En.size() << std::endl;

    Timer t2;    
    t2.Start();

    int total = 0;
#pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
    for (NodeID u=0; u < SISAg.num_nodes(); u++) {
        for (NodeID v : SISAg.neigh(u)) {
            if (v > u)
                break;
            auto it = SISAg.neigh(u).begin();
            for (NodeID w : SISAg.neigh(v)) {
                if (w > v)
                    break;
                while (*it < w)
                    it++;
                if (w == *it)
                    total++;
            }
        }
    }
    total = (int) ((1.*total)*(1./p)*(1./p)*(1./p));


    t2.Stop();
    PrintTime("Counting time", t2.Seconds());
    tc_time = t2.Seconds();

    auto threshold = 0;
    auto k = 0;
    auto m = 0;
    auto approx_str_size = 0;
    auto initial_csr_size = g.getSize() / (1024.0 * 1024.0); //MB

    // RRR - this means that a given line is dedicated to the runtime results
    // the columns are as follows:
    // RRR [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [preprocessing-time] [tc-time] [total-runtime] [approximated TC count]

    std::cout << "RRR TC DOULION TC_DOULION " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " "  << m << " " << pp_time << " " << tc_time << " " << pp_time + tc_time <<  " " << total << std::endl;

    // SSS - this means that a given line is dedicated to the size results
    // the columns are as follows:
    // SSS [Problem] [approximation-scheme] [baseline (problem + approx-scheme)] [graph-name] [thread count] [number of vertices] [number of edges] [threshold (KMV parameter)] [k (another KMV parameter-number of hash functions] [size of BF structures] [size of the original standard CSR (total)] [total size of both] [approximated TC count]

    std::cout << "SSS TC DOULION TC_DOULION " << graphName << " " << threads << " " << g.num_nodes() << " " << g.num_edges() << " " << threshold << " " << k << " " << m << " " << approx_str_size << " " << initial_csr_size << " " << approx_str_size + initial_csr_size << " " << total << std::endl;



    return total;


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


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
    cout << total_triangles << " triangles" << endl;
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
        return TCDoulion(g, cli.p(), cli.getGraphBasename(), cli.getThreadNum());
    };
    BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);



    return 0;
}

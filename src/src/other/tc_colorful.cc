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

int TCColourful(const Graph &g, float  p_)
{


    SISA_Graph SISAg(g);

    std::random_device rd;
    std::mt19937 eng(22);//rd());

    SSSet E = SISAg.GetEdges();
    int64_t n = SISAg.num_nodes();

    int N = 1/p_;

    std::uniform_int_distribution<> distr(1,N);

    vector<int> colors(n);
#pragma omp parallel for
    for (NodeID u=0; u < SISAg.num_nodes(); u++)
        colors[u] = distr(eng);


#pragma omp parallel for
    for(auto e: E)
        if (E2U(e,n) < E2V(e,n) and colors[E2U(e,n)] != colors[E2V(e,n)])
        {
                SISAg.RemoveEdge(E2U(e,n) , E2V(e,n));
                SISAg.RemoveEdge(E2V(e,n) , E2U(e,n));
        }

//    SSSet E_prime = SISAg.GetEdges();
//        for (auto e: E_prime)
//            if(E2U(e,n) < E2V(e,n))
//                cout<<E2U(e,n)<<" - "<<E2V(e,n)<<endl;
//
//    for (NodeID u=0; u < SISAg.num_nodes(); u++) {
//        cout<<"Node: "<<u<<": ";
//        for (NodeID v : SISAg.neigh(u)) {
//            cout<<v<<" ";
//        }
//        cout<<endl;
//    }

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


    return total*N*N;


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
        return TCColourful(g, cli.p());
    };
    BenchmarkKernel(cli, g, ct, PrintTriangleStats, TCVerifier);



    return 0;
}

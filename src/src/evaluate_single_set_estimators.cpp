#include <iostream>
#include <fstream>
#include <string>
#include <math.h> 
#include <algorithm>
#include "sets.hpp"
#include "MurmurHash3.h"

// Read a given graph. Return both, its adjacency list and a list of its edges
std::tuple<std::vector<std::vector<int>>,std::vector<std::tuple<int,int>>> read_graph(std::string name){
	std::string path = "../graphs/" + name + ".el";
    std::fstream graphfile(path, std::ios_base::in);
	// a and b are the vertex-ids of the two vertices adjacent to an edge
    int a, b = 0;
	// We find the number of vertices n by looking for the maximum vertex id and adding 1
	int n = 0;
	// Read all edges
	std::vector<std::tuple<int,int>> edges;
    while (graphfile >> a >> b)
    {
		edges.push_back(std::make_tuple(a,b));
		n = n > a ? n : a;
		n = n > b ? n : b;
    }
	n += 1;
	// Create adjacency list
    std::vector<std::vector<int>> adj_list(n, std::vector<int>());
	for(std::tuple<int,int> edge : edges){
		a = std::get<0>(edge);
		b = std::get<1>(edge);
		adj_list[a].push_back(b);
		adj_list[b].push_back(a);
	}
	// Sort neighborhoods, this is needed for the computation of the correct intersection
	// set intersection cardinality.
	for(int i = 0; i < adj_list.size(); i++){
		std::sort(adj_list[i].begin(), adj_list[i].end());
	}
	// Return adjacency list and edge list
	return std::make_tuple(adj_list, edges);
}

int estimator_1(int b, int Bx, int Bx1){
    if(Bx == 0)
        return 0;
    float est = (float)Bx1 / (float)b;
    return (int)(est + 0.5);
}

int estimator_2(int b, int Bx, int Bx1){
    if(Bx == 0)
        return 0;
    float est = -1 * ((float)Bx / (float)b) * log(1 - ((float)Bx1 / float(Bx)));
    return static_cast<unsigned int>(est + 0.5);
}


// Evaluate a given graph for all combinations of parameters b in all_b and s in all_s
int evaluate_estimators_on_graph(std::string name, std::vector<int> all_b, std::vector<int> all_s){
	// Read graph
	std::cout << "Reading Graph" << std::endl;
	std::tuple<std::vector<std::vector<int>>,std::vector<std::tuple<int,int>>> graph = read_graph(name);
	std::vector<std::vector<int>> adj_list= std::get<0>(graph);
	std::vector<std::tuple<int,int>> edges = std::get<1>(graph);
	// Compute graph size in bits
	int n = adj_list.size();
	int m = 0;
	for(int i = 0; i < n; i++){
		m += adj_list[i].size();
	}
	m /= 2;
	float graph_size = n * 64 + 2 * m * 64; // in bits

	// Evaluate estimators for different values of b and s
	for (int b : all_b){
		for (int s : all_s){
			std::cout << "Computing Estimates with b = " << b << " and s = " << s << std::endl;
			// Compute number of bits in the bloom filter and number of hash functions in 1-hash / k-hash
			float s_f = (100.0 + (float)s) / 100.0;							// Storage (factor)
			int bf_size = ((int)(graph_size * (s_f - 1.0) / n / 8.0))*8;	// Bits in bloom filter
			int k = (int)(graph_size * (s_f - 1.0) / n / 64);				// Hashes in 1-hash / k-hash
			// Construct Seeds 
			struct timeval seed;
			gettimeofday(&seed, NULL);
			srand(seed.tv_usec);
			std::vector<uint32_t> seeds;
			for (int i = 0; i<b; i++) {
				seeds.push_back(rand());
			}

			// Open file to store results
			std::string type = all_b.size() == 1 ? "vm" : "vb";
			std::string out_path = "../../patrick/plots/results_estimators/single_set_" + type + "_" + name + "_" + std::to_string(b) + "_" + std::to_string(s) + ".csv";
			std::ofstream outfile;
			outfile.open(out_path);

			std::cout << "Computing Results" << std::endl;
			// Evaluate both estimators for all neighborhoods
			for(int i = 0; i < adj_list.size(); i++){
				std::vector<int> elements;
				for(int j = 0; j < adj_list[i].size(); j++){
					elements.push_back(adj_list[i][j]);
				}

				int* start = &*elements.begin();
				int* end = &*elements.end();
				int size = adj_list.size();
				int bb = (all_b.size() == 1) ? 1 : b;
				BloomSet<int> * set_1 = new BloomSet<int>(start, end, bf_size, bb, seeds, size);
				BloomSet<int> * set_b = new BloomSet<int>(start, end, bf_size, b, seeds, size);

				int Bx1_1 = set_1->bitarray.count();
				int Bx1_b = set_b->bitarray.count();
				int Bx = bf_size;

                int correct = adj_list[i].size();
                int est_1 = estimator_1(bb, Bx, Bx1_1);
                int est_2 = estimator_2(b, Bx, Bx1_b);

				outfile << correct << ", " << est_1 << "," << est_2 << std::endl;
			}
			outfile.close();
		}
	}
	return 0;
}

int main(){
	std::string config_file= "../../patrick/plots/config/graphs_single_set.csv";
    std::fstream cfgfile(config_file, std::ios_base::in);

	std::string graph_name;
    while (cfgfile >> graph_name){
		std::cout << graph_name << std::endl;
	}

	std::vector<int> all_b;
	std::vector<int> all_s;

	// Vary b
	all_b = {1,2,3,4,5,6,7,8};
	all_s = {33};
	evaluate_estimators_on_graph(graph_name, all_b, all_s);

	// Vary memory
	all_b = {3};
	all_s = {10,20,25,33,50};
	evaluate_estimators_on_graph(graph_name, all_b, all_s);

    return 0;
}

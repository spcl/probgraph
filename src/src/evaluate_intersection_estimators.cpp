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

	// Compute the correct values for set intersections
	std::cout << "Computing Correct Values" << std::endl;
	std::map<std::tuple<int,int>,int> correct_values;
	for (std::tuple<int,int> edge : edges){
		std::vector<int> A = adj_list[std::get<0>(edge)];
		std::vector<int> B = adj_list[std::get<1>(edge)];
		std::vector<int> C;
		std::set_intersection(A.begin(),A.end(),B.begin(),B.end(),back_inserter(C));
		correct_values[edge] = C.size();
	}

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
			std::string out_path = "../../intersection_estimator_results/intersection_" + name + "_" + std::to_string(b) + "_" + std::to_string(s) + ".csv";
			std::ofstream outfile;
			outfile.open(out_path);

			// Construct approximate set representations for each neighborhood
			std::vector<BloomSet<int>*>	bloomsets;
			std::vector<OneHashSet<int>*> onehashsets;
			std::vector<KMinHashSet<int>*> khashsets;
			for(int i = 0; i < n; i++){
				int* start = &(*(adj_list[i].begin()));
				int* end = &(*(adj_list[i].end()));
				bloomsets.push_back(new BloomSet<int>(start, end, bf_size, b, seeds, graph_size));
				onehashsets.push_back(new OneHashSet<int>(start, end, k, seeds[0]));
				khashsets.push_back(new KMinHashSet<int>(start, end, k, seeds[0]));
			}
			std::cout << "Computing Results" << std::endl;
			// Evaluate all 5 estimators for all edges
			for (std::tuple<int,int> edge : edges){
				int aa = std::get<0>(edge);					// vertex id of first vertex
				int bb = std::get<1>(edge);					// vertex id of second vertex
				BloomSet<int>* bs1 = bloomsets[aa];
				BloomSet<int>* bs2 = bloomsets[bb];
				OneHashSet<int>* ohs1 = onehashsets[aa];
				OneHashSet<int>* ohs2 = onehashsets[bb];
				KMinHashSet<int>* khs1 = khashsets[aa];
				KMinHashSet<int>* khs2 = khashsets[bb];
				std::vector<int> A = adj_list[aa];
				std::vector<int> B = adj_list[bb];

				int correct = correct_values[edge];
				int est_1 = bs1->intersect_count(*bs2);
				//int est_2 = A.size() + B.size() - bs1->union_count(*bs2);
				int est_3 = (bs1->bitarray & bs2->bitarray).count() / b;
				int est_4 = ohs1->intersect_count(*ohs2);
				int est_5 = khs1->intersect_count(*khs2);
				//outfile << correct << ", " << est_1 << "," << est_2 << ", " << est_3 <<  "," << est_4 << ", " << est_5 << std::endl;
				outfile << correct << ", " << est_1 << "," << est_3 <<  "," << est_4 << ", " << est_5 << std::endl;
			}
			outfile.close();
		}
	}
	return 0;
}

int main(){
	std::string config_file= "../../intersection_estimator_config/graphs_intersection.csv";
    std::fstream cfgfile(config_file, std::ios_base::in);

	std::vector<int> all_b = {1,4};
	std::vector<int> all_s = {33};

	std::string graph_name;
    while (cfgfile >> graph_name){
		std::cout << graph_name << std::endl;
		evaluate_estimators_on_graph(graph_name, all_b, all_s);
	}
    return 0;
}

#include <string>
#include <iostream>
#include <fstream>
#include "fm.h"

using namespace hypergraph;

double parse_epsilon(std::string file_name) {
    std::fstream graph_file(file_name);
    assert(graph_file.is_open());

    double epsilon;
    graph_file >> epsilon;
    return epsilon;
}

HyperGraph<int, int> parse_hypergraph(std::string file_name) {
    std::fstream graph_file(file_name);
    assert(graph_file.is_open());

    double epsilon;
    std::vector<std::vector<IndexType>> edge2node_map;

    // parse graph file
    graph_file >> epsilon;
    IndexType max_node_id = 0;
    std::string token;
    while (graph_file >> token)
    {
        assert (token == "NET");
        graph_file >> token;
        std::cout << "Parse Net: " << token << std::endl;
        std::vector<IndexType> edge;
        graph_file >> token;
        while (token != ";") {
            assert(token[0] == 'c');
            
            IndexType cell_sufix = std::stoi(token.substr(1));
            IndexType node_id = cell_sufix - 1;
            max_node_id = std::max(max_node_id, node_id);
            edge.push_back(node_id);

            graph_file >> token;
        }

        edge2node_map.push_back(edge);
    }

    // build hypergraph

    HyperGraph<int, int> hypergraph;
    for (size_t i = 0; i <= max_node_id; i++) {
        hypergraph.add_node(1);
    }

    for (const std::vector<IndexType>& edge : edge2node_map) {
        hypergraph.add_edge(1, edge);
    }

    hypergraph.build_node2edge_map();

    return hypergraph;
}

int main() {
    std::string file_path = "./benchmarks/partition/simple.dat";
    HyperGraph<int, int> hypergraph = parse_hypergraph(file_path);
    double epsilon = parse_epsilon(file_path);
    fm_partitioner::FMPartitioner<int, int> partitioner(hypergraph, epsilon);
    
    partitioner();
    return 0;
}
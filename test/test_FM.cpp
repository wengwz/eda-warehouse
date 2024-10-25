#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
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
        //std::cout << "Parse Net: " << token << std::endl;
        std::set<IndexType> edge_nodes_set;
        graph_file >> token;
        while (token != ";") {
            assert(token[0] == 'c');
            
            IndexType cell_sufix = std::stoi(token.substr(1));
            IndexType node_id = cell_sufix - 1;
            max_node_id = std::max(max_node_id, node_id);
            edge_nodes_set.insert(node_id);

            graph_file >> token;
        }

        std::vector<IndexType> edge_nodes(edge_nodes_set.begin(), edge_nodes_set.end());
        edge2node_map.push_back(edge_nodes);
    }

    // build hypergraph

    HyperGraph<int, int> hypergraph;
    for (size_t i = 0; i <= max_node_id; i++) {
        hypergraph.add_node(1);
    }

    for (const std::vector<IndexType>& edge : edge2node_map) {
        if (edge.size() < 2) continue;
        hypergraph.add_edge(1, edge);
    }

    hypergraph.build_node2edge_map();

    return hypergraph;
}

void write_hypergraph_hmetis(std::string file_name, const HyperGraph<int, int>& graph) {
    std::ofstream hmetis_graph_file(file_name);
    assert(hmetis_graph_file.is_open());
    hmetis_graph_file << graph.get_edges_num() << " " << graph.get_nodes_num() << std::endl;

    for (IndexType edgeId = 0; edgeId < graph.get_edges_num(); edgeId++) {
        const std::vector<IndexType>& edge_nodes = graph.get_nodes_of_edge(edgeId);
        for (int i = 0; i < edge_nodes.size(); i++) {
            if (i > 0) {
                hmetis_graph_file << " ";
            }
            hmetis_graph_file << edge_nodes[i] + 1;
        }      
        hmetis_graph_file << std::endl;
    }

    hmetis_graph_file.close();
}

int main() {
    std::string file_path = "./benchmarks/partition/input_0.dat";
    HyperGraph<int, int> hypergraph = parse_hypergraph(file_path);
    write_hypergraph_hmetis("./benchmarks/partition/input_0.hgr", hypergraph);
    double epsilon = parse_epsilon(file_path);
    partition::FMPartitioner<int, int> partitioner(hypergraph, epsilon);
    
    partitioner();
    return 0;
}
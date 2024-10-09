#ifndef _HYPER_GRAPH_H_
#define _HYPER_GRAPH_H_

#include <vector>
#include <cstdint>
#include <cassert>
#include <algorithm>

namespace hypergraph {

using IndexType = int;


template<typename WeightType>
struct Node
{
    std::vector<IndexType> edge_ids;
    WeightType weight;
};

template<typename WeightType>
struct HyperEdge
{
    WeightType weight;
    std::vector<IndexType> node_ids;

};

template<typename NodeWeightType, typename EdgeWeightType>
class HyperGraph {
    using GraphNode = Node<NodeWeightType>;
    using GraphEdge = HyperEdge<EdgeWeightType>;

public:
    HyperGraph() : num_nodes(0), num_edges(0) {};

    IndexType add_node(const NodeWeightType& weight) {
        GraphNode node;
        node.weight = weight;
        nodes.push_back(node);
        return num_nodes++;
    }

    IndexType add_edge(const EdgeWeightType& weight, const std::vector<IndexType>& node_ids) {
        assert (node_ids.size() >= 0);
        IndexType max_node_id = *(std::max_element(node_ids.begin(), node_ids.end()));
        if (max_node_id >= num_nodes) {
            throw std::runtime_error("Node id out of range");
        }

        GraphEdge edge;
        edge.weight = weight;
        edge.node_ids = node_ids;
        edges.push_back(edge);
        return num_edges++;
    }

    void build_node2edge_map() {
        assert(num_nodes > 0 && num_edges > 0);

        for (GraphNode& node : nodes) {
            node.edge_ids.clear();
        }

        for (IndexType i = 0; i < num_edges; i++) {
            for (IndexType nodeId : edges[i].node_ids) {
                nodes[nodeId].edge_ids.push_back(i);
            }
        }

        // check nodes
        for (IndexType i = 0; i < num_nodes; i++) {
            if (nodes[i].edge_ids.empty()) {
                throw std::runtime_error("Node-" + std::to_string(i) + " has no incident edges");
            }
        }
    }

    // getter
    IndexType get_nodes_num() const {
        return num_nodes;
    }

    IndexType get_edges_num() const {
        return num_edges;
    }

    IndexType get_degree_of_edge(IndexType edgeId) const {
        return edges[edgeId].node_ids.size();
    }

    const std::vector<IndexType>& get_nodes_of_edge(IndexType edgeId) const {
        return edges[edgeId].node_ids;
    }

    const std::vector<IndexType>& get_edges_of_node(IndexType nodeId) const {
        return nodes[nodeId].edge_ids;
    }

    NodeWeightType get_weight_of_node(IndexType nodeId) const {
        return nodes[nodeId].weight;
    }

    EdgeWeightType get_weight_of_edge(IndexType edgeId) const {
        return edges[edgeId].weight;
    }

    NodeWeightType get_total_node_weight() const {
        NodeWeightType total_weight = 0;
        for (const GraphNode& node : nodes) {
            total_weight += node.weight;
        }
        return total_weight;
    }

    EdgeWeightType get_total_edge_weight() const {
        EdgeWeightType total_weight = 0;
        for (const GraphEdge& edge : edges) {
            total_weight += edge.weight;
        }
        return total_weight;
    }

private:
    IndexType num_nodes;
    IndexType num_edges;
    std::vector<GraphNode> nodes;
    std::vector<GraphEdge> edges;
};

}

#endif
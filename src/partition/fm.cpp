#include <iostream>
#include <sstream>
#include <cmath>
#include "fm.h"

namespace fm_partitioner {

template<typename NodeWeightType, typename EdgeWeightType>
FMPartitioner<NodeWeightType, EdgeWeightType>::FMPartitioner(
    const std::vector<NodeWeightType>& node_weights,
    const std::vector<EdgeWeightType>& edge_weights,
    const std::vector<std::vector<IndexType>>& edge2node, 
    double epsilon
) {
    this->epsilon = epsilon;
    for (const auto& weight : node_weights) {
        hypergraph.add_node(weight);
    }
    assert (edge_weights.size() == edge2node.size());
    for (size_t i = 0; i < edge_weights.size(); i++) {
        hypergraph.add_edge(edge_weights[i], edge2node[i]);
    }
    hypergraph.build_node2edge_map();
}

template<typename NodeWeightType, typename EdgeWeightType>
void FMPartitioner<NodeWeightType, EdgeWeightType>::operator()() {
    std::stringstream ss;
    ss << "Start 2-way partition with";
    ss << " epsilon=" << epsilon;
    ss << " node_num=" << hypergraph.get_nodes_num();
    ss << " edge_num=" << hypergraph.get_edges_num();

    std::cout << ss.str() << std::endl;

    BlkSizeType max_blk_size = get_max_blk_size();
    BlkSizeType min_blk_size = get_min_blk_size();

    init_partition();
    int passId = 0;
    EdgeWeightType cut_size = get_cutsize();
    std::array<BlkSizeType, 2> blk_sizes = {get_size_of_blk(0), get_size_of_blk(1)};

    std::cout << "Initial partition: " << get_partition_info(cut_size, blk_sizes) << std::endl;

    while (true)
    {
        MoveGainType pass_gain = 0;
        MoveGainType pass_max_gain = 0;
        IndexType pass_max_gain_step = -1;
        std::vector<IndexType> pass_node_move_steps;
        
        std::vector<bool> is_node_moved(hypergraph.get_nodes_num(), false);
        std::array<BlkSizeType, 2> tmp_blk_sizes = blk_sizes;
        std::vector<IndexType> tmp_node_blk_ids = node_blk_ids;

        while (true)
        {// tentatively move nodes to find pass_max_gain and pass_max_gain_step
            MoveGainType node_max_gain = 0;
            IndexType max_gain_nodeId = -1;

            for (IndexType nodeId = 0; nodeId < hypergraph.get_nodes_num(); ++nodeId) {
                // find node of max gain if moved
                //std::cout << "check node-" << nodeId << std::endl;
                if (is_node_moved[nodeId]) {
                    continue;
                }

                // check block size constraint
                NodeWeightType node_weight = hypergraph.get_weight_of_node(nodeId);
                IndexType from_blk_id = tmp_node_blk_ids[nodeId];
                IndexType to_blk_id = get_opposite_blk_id(from_blk_id);
                BlkSizeType from_blk_size = tmp_blk_sizes[from_blk_id] - node_weight;
                BlkSizeType to_blk_size = tmp_blk_sizes[to_blk_id] + node_weight;
                if (from_blk_size < min_blk_size || to_blk_size > max_blk_size) {
                    continue;
                }

                MoveGainType node_gain = get_gain_of_node(nodeId, tmp_node_blk_ids);
                if (node_gain > node_max_gain || max_gain_nodeId == -1) {
                    node_max_gain = node_gain;
                    max_gain_nodeId = nodeId;
                }
            }

            if (max_gain_nodeId == -1) {
                break; // finish this pass if unable to find nodes to move
            }

            // tentatively move node and record the nodeId
            is_node_moved[max_gain_nodeId] = true;
            NodeWeightType node_weight = hypergraph.get_weight_of_node(max_gain_nodeId);
            tmp_blk_sizes[tmp_node_blk_ids[max_gain_nodeId]] -= node_weight;
            IndexType to_blk_id = get_opposite_blk_id(tmp_node_blk_ids[max_gain_nodeId]);
            tmp_blk_sizes[to_blk_id] += node_weight;
            tmp_node_blk_ids[max_gain_nodeId] = to_blk_id;


            pass_gain += node_max_gain;
            pass_node_move_steps.push_back(max_gain_nodeId);

            if (pass_gain > pass_max_gain || pass_max_gain_step == -1) {
                pass_max_gain = pass_gain;
                pass_max_gain_step = pass_node_move_steps.size();
            }
            std::cout << "find max-gain node=" << max_gain_nodeId << " gain=" << node_max_gain << std::endl;
        }

        if (pass_max_gain == 0) {
            break; // finish iterative node moving if no more gain
        }

        // move nodes within the pass_max_gain_step
        for (IndexType i = 0; i < pass_max_gain_step; i++) {
            IndexType nodeId = pass_node_move_steps[i];
            NodeWeightType node_weight = hypergraph.get_weight_of_node(nodeId);

            IndexType from_blk_id = node_blk_ids[nodeId];
            blk_sizes[from_blk_id] -= node_weight;
            IndexType to_blk_id = get_opposite_blk_id(from_blk_id);
            blk_sizes[to_blk_id] += node_weight;

            node_blk_ids[nodeId] = tmp_node_blk_ids[nodeId];
        }
        cut_size -= pass_max_gain;

        std::cout << "Pass-" << passId << ": " << get_partition_info(cut_size, blk_sizes) << std::endl;
        passId += 1;
    }

}

// getter functions
template<typename NodeWeightType, typename EdgeWeightType>
EdgeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_cutsize() const {

    if (node_blk_ids.size() != hypergraph.get_nodes_num()) {
        std::runtime_error("Unable to get cut size due to invalid partition results");
    }

    EdgeWeightType cut_size = 0;

    for (IndexType edgeId = 0; edgeId < hypergraph.get_edges_num(); ++edgeId) {
        const std::vector<IndexType> edge_nodes = hypergraph.get_nodes_of_edge(edgeId);
        EdgeWeightType edge_weight = hypergraph.get_weight_of_edge(edgeId);
        
        uint32_t block_id = get_blk_id_of_node(edge_nodes[0]);
        for (const IndexType& nodeId : edge_nodes) {
            if (block_id != get_blk_id_of_node(nodeId)) {
                cut_size += edge_weight;
                break;
            }
        }
    }

    return cut_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_size_of_blk(IndexType blkId) const {
    if (blkId >= 2) {
        std::runtime_error("Invalid block id for 2-way partition");
    }

    if (node_blk_ids.size() != hypergraph.get_nodes_num()) {
        std::runtime_error("Unable to get block size due to invalid partition results");
    }

    BlkSizeType blk_size = 0;
    for (IndexType nodeId = 0; nodeId < hypergraph.get_nodes_num(); ++nodeId) {
        if (node_blk_ids[nodeId] == blkId) {
            blk_size += hypergraph.get_weight_of_node(nodeId);
        }
    }

    return blk_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
std::vector<IndexType> FMPartitioner<NodeWeightType, EdgeWeightType>::get_nodes_of_blk(IndexType blkId) const {
    if (blkId >= 2) {
        std::runtime_error("Invalid block id for 2-way partition");
    }

    if (node_blk_ids.size() != hypergraph.get_nodes_num()) {
        std::runtime_error("Unable to get nodes of block due to invalid partition results");
    }

    std::vector<IndexType> node_ids;
    
    for (IndexType nodeId = 0; nodeId < hypergraph.get_nodes_num(); ++nodeId) {
        if (node_blk_ids[nodeId] == blkId) {
            node_ids.push_back(nodeId);
        }
    }

    return node_ids;
}

template<typename NodeWeightType, typename EdgeWeightType>
IndexType FMPartitioner<NodeWeightType, EdgeWeightType>::get_blk_id_of_node(IndexType node) const {
    if (node >= hypergraph.get_nodes_num()) {
        std::runtime_error("Invalid node id");
    }

    return node_blk_ids[node];
}

template<typename NodeWeightType, typename EdgeWeightType>
const std::vector<IndexType>& FMPartitioner<NodeWeightType, EdgeWeightType>::get_blk_ids() const {
    return node_blk_ids;
}


// private helper functions
template<typename NodeWeightType, typename EdgeWeightType>
void FMPartitioner<NodeWeightType, EdgeWeightType>::init_partition() {
    // TODO:
    node_blk_ids.clear();
    for (IndexType i = 0; i < hypergraph.get_nodes_num(); ++i) {
        uint32_t part_res = (i % 2 == 0) ? 0 : 1;
        node_blk_ids.push_back(part_res);
    }
}

template<typename NodeWeightType, typename EdgeWeightType>
EdgeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_gain_of_node(IndexType nodeId, const std::vector<IndexType>& tmp_node_blk_ids) const {
    MoveGainType node_gain = 0;
    IndexType node_blk_id = tmp_node_blk_ids[nodeId];
    for (const IndexType& edgeId : hypergraph.get_edges_of_node(nodeId)) {
        IndexType num_same_blk_nodes = 0;
        IndexType edge_degree = hypergraph.get_degree_of_edge(edgeId);
        
        for (IndexType edge_nodeId : hypergraph.get_nodes_of_edge(edgeId)) {
            if (tmp_node_blk_ids[edge_nodeId] == node_blk_id) {
                num_same_blk_nodes += 1;
            }
        }

        if (num_same_blk_nodes == 1) {
            node_gain += hypergraph.get_weight_of_edge(edgeId);
        } else if (num_same_blk_nodes == edge_degree) {
            node_gain -= hypergraph.get_weight_of_edge(edgeId);
        }
    }
    return node_gain;
}

template<typename NodeWeightType, typename EdgeWeightType>
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_size_of_blk(IndexType blkId, const std::vector<IndexType>& tmp_node_blk_ids) const {
    BlkSizeType blk_size = 0;
    for (IndexType nodeId = 0; nodeId < hypergraph.get_nodes_num(); ++nodeId) {
        if (tmp_node_blk_ids[nodeId] == blkId) {
            blk_size += hypergraph.get_weight_of_node(nodeId);
        }
    }
    return blk_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_max_blk_size() const {
    NodeWeightType total_node_weight = hypergraph.get_total_node_weight();
    double half_node_weight = (double) total_node_weight / 2.0;
    return std::floor(half_node_weight * (1 + epsilon));
}

template<typename NodeWeightType, typename EdgeWeightType>
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_min_blk_size() const {
    BlkSizeType total_node_weight = hypergraph.get_total_node_weight();
    double half_node_weight = (double) total_node_weight / 2.0;
    return std::ceil(half_node_weight * (1 - epsilon));
}

template<typename NodeWeightType, typename EdgeWeightType>
std::string FMPartitioner<NodeWeightType, EdgeWeightType>::get_partition_info(EdgeWeightType cut_size, const std::array<BlkSizeType, 2>& blk_sizes) {
    std::stringstream ss;
    ss << "cut_size=" << cut_size;
    ss << " blk0=" << blk_sizes[0];
    ss << " blk1=" << blk_sizes[1];
    return ss.str();
}

} // namespace fm_partitioner
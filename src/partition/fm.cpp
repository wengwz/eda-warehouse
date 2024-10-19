#include <iostream>
#include <sstream>
#include <cmath>
#include <stack>
#include "fm.h"

namespace partition {

template<typename NodeWeightType, typename EdgeWeightType>
FMPartitioner<NodeWeightType, EdgeWeightType>::FMPartitioner(
    const std::vector<NodeWeightType>& node_weights,
    const std::vector<EdgeWeightType>& edge_weights,
    const std::vector<std::vector<IndexType>>& edge2node, 
    double epsilon
):epsilon(epsilon), gain_bucket_array(node_weights.size()) {

    // build hypergraph
    for (const auto& weight : node_weights) {
        hypergraph.add_node(weight);
    }
    assert (edge_weights.size() == edge2node.size());
    for (size_t i = 0; i < edge_weights.size(); i++) {
        if (edge2node[i].size() < 2) continue;
        hypergraph.add_edge(edge_weights[i], edge2node[i]);
    }
    hypergraph.build_node2edge_map();

    // init partition
    init_partition();
}

template<typename NodeWeightType, typename EdgeWeightType>
void FMPartitioner<NodeWeightType, EdgeWeightType>::run() {
    std::stringstream ss;
    ss << "Start 2-way partition with";
    ss << " epsilon=" << epsilon;
    ss << " node_num=" << hypergraph.get_nodes_num();
    ss << " edge_num=" << hypergraph.get_edges_num();

    std::cout << ss.str() << std::endl;

    int passId = 0;
    std::cout << "Initial partition: " << get_partition_info() << std::endl;

    while (true)
    { // one pass
        // record pass gain during trial moves
        MoveGainType pass_gain = 0;
        MoveGainType pass_max_gain = 0;
        IndexType pass_max_gain_step = -1;
        std::vector<IndexType> pass_node_move_steps;
        
        // partition states
        std::array<BlkSizeType, 2> tmp_blk_sizes = blk_sizes;
        std::vector<IndexType> tmp_node_blk_ids = node_blk_ids;
        BucketArrayType tmp_gain_bucket_array(hypergraph.get_nodes_num());
        // init tmp_gain_bucket_array (TODO: to be removed)
        for (IndexType nodeId = 0; nodeId < hypergraph.get_nodes_num(); nodeId++) {
            MoveGainType node_gain = get_gain_of_node(nodeId);
            tmp_gain_bucket_array.insert(nodeId, node_gain);
        }

        while (true)
        {// tentatively move nodes to find pass_max_gain and pass_max_gain_step

            //search feasible node with max gain
            std::stack<std::pair<IndexType, MoveGainType>> illegal_nodes;
            IndexType nodeId;
            MoveGainType node_gain;
            while (true)
            {
                nodeId = tmp_gain_bucket_array.get_top_id();
                if (nodeId == -1) { // no more nodes in the bucket array
                    break;
                }
                node_gain = tmp_gain_bucket_array.get_value_of(nodeId);
                tmp_gain_bucket_array.pop();

                if (is_move_legal(nodeId, tmp_node_blk_ids[nodeId], tmp_blk_sizes)) {
                    break;
                }
                illegal_nodes.push(std::make_pair(nodeId, node_gain));
            }

            if (nodeId == -1) {
                break; // finish pass if no more feasible nodes in tmp_gain_bucket_array
            }

            while (!illegal_nodes.empty())
            {
                IndexType illegal_nodeId = illegal_nodes.top().first;
                MoveGainType illegal_node_gain = illegal_nodes.top().second;
                illegal_nodes.pop();
                tmp_gain_bucket_array.insert(illegal_nodeId, illegal_node_gain);
            }

            move_node(nodeId, tmp_node_blk_ids, tmp_blk_sizes, tmp_gain_bucket_array);
            
            pass_gain += node_gain;
            pass_node_move_steps.push_back(nodeId);
            if (pass_gain > pass_max_gain || pass_max_gain_step == -1) {
                pass_max_gain = pass_gain;
                pass_max_gain_step = pass_node_move_steps.size();
            }

            //std::cout << "The size of tentatively moved nodes: " << pass_node_move_steps.size() << std::endl;
            //std::cout << "Tentativly move node-" << nodeId << " with gain=" << node_gain << std::endl;
        }

        if (pass_max_gain <= 0) {
            break; // finish iterative node moving if no more gain
        }

        // move nodes within the pass_max_gain_step
        for (IndexType step = 0; step < pass_max_gain_step; step++) {
            IndexType nodeId = pass_node_move_steps[step];
            cut_size -= get_gain_of_node(nodeId);
            // std::cout << "Move node-" << nodeId << " with gain=" << get_gain_of_node(nodeId) << " cut_size=" << cut_size << std::endl;
            move_node(nodeId, node_blk_ids, blk_sizes, gain_bucket_array);

            // if (nodeId == 144565 || nodeId == 144579 || nodeId == 144587) {
            //     std::vector<MoveGainType> nodes_gain = get_gain_of_nodes(node_blk_ids);
            //     for (IndexType i = 0; i < hypergraph.get_nodes_num(); i++) {
            //         MoveGainType node_gain_ref = nodes_gain[i];
            //         MoveGainType node_gain_incr = get_gain_of_node(i);
            //         if (node_gain_incr != node_gain_ref) {
            //             std::cout << "Node-" << i << " gain=" << node_gain_ref << " gain_incr=" << node_gain_incr << std::endl;
            //         }
            //     }
            // }

            // EdgeWeightType actual_cut_size = get_cutsize(node_blk_ids);
            // if (actual_cut_size != cut_size) {
            //     std::cout << "Actual cut size: " << actual_cut_size << std::endl;
            //     assert(actual_cut_size == cut_size);
            // }
        }


        std::cout << "Pass-" << passId << ": " << get_partition_info() << std::endl;
        //std::cout << "Actual cut size: " << get_cutsize(node_blk_ids) << std::endl;
        passId += 1;
    }

}

// getter functions
template<typename NodeWeightType, typename EdgeWeightType>
EdgeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_cutsize() const {
    return cut_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_size_of_blk(IndexType blkId) const {
    if (blkId >= 2) {
        std::runtime_error("Invalid block id for 2-way partition");
    }

    return blk_sizes[blkId];
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
const std::vector<IndexType>& FMPartitioner<NodeWeightType, EdgeWeightType>::get_partition_result() const {
    return node_blk_ids;
}


// private helper functions
template<typename NodeWeightType, typename EdgeWeightType>
void FMPartitioner<NodeWeightType, EdgeWeightType>::init_partition() {
    // init partition results
    // TODO:
    node_blk_ids.clear();
    for (IndexType i = 0; i < hypergraph.get_nodes_num(); ++i) {
        // assign initial partition result according to nodeId
        uint32_t part_res = (i % 2 == 0) ? 0 : 1;
        node_blk_ids.push_back(part_res);
    }
    cut_size = get_cutsize(node_blk_ids);
    blk_sizes[0] = get_size_of_blk(0, node_blk_ids);
    blk_sizes[1] = get_size_of_blk(1, node_blk_ids);

    // init gain bucket array
    std::vector<MoveGainType> node_gains = get_gain_of_nodes(node_blk_ids);
    for (IndexType nodeId = 0; nodeId < hypergraph.get_nodes_num(); nodeId++) {
        MoveGainType node_gain = node_gains[nodeId];
        gain_bucket_array.insert(nodeId, node_gain);
    }

    //
    blk_max_size = get_blk_max_size();
    blk_min_size = get_blk_min_size();
}

template<typename NodeWeightType, typename EdgeWeightType>
EdgeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_cutsize(const std::vector<IndexType> tmp_node_blk_ids) const {

    assert(tmp_node_blk_ids.size() == hypergraph.get_nodes_num());

    EdgeWeightType tmp_cut_size = 0;

    for (IndexType edgeId = 0; edgeId < hypergraph.get_edges_num(); ++edgeId) {
        const std::vector<IndexType>& edge_nodes = hypergraph.get_nodes_of_edge(edgeId);
        EdgeWeightType edge_weight = hypergraph.get_weight_of_edge(edgeId);
        
        uint32_t block_id = tmp_node_blk_ids[edge_nodes[0]];
        for (const IndexType& nodeId : edge_nodes) {
            if (block_id != tmp_node_blk_ids[nodeId]) {
                tmp_cut_size += edge_weight;
                break;
            }
        }
    }

    return tmp_cut_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
std::vector<EdgeWeightType> FMPartitioner<NodeWeightType, EdgeWeightType>::get_gain_of_nodes(const std::vector<IndexType>& tmp_node_blk_ids) const {
    std::vector<MoveGainType> node_gains(hypergraph.get_nodes_num(), 0);

    for (IndexType edgeId = 0; edgeId < hypergraph.get_edges_num(); edgeId++) {
        EdgeWeightType edge_weight = hypergraph.get_weight_of_edge(edgeId);

        std::array<std::vector<IndexType>, 2> blk2node;
        for (IndexType nodeId : hypergraph.get_nodes_of_edge(edgeId)) {
            IndexType blkId = tmp_node_blk_ids[nodeId];

            blk2node[blkId].push_back(nodeId);
        }

        if (blk2node[0].size() == 1 || blk2node[1].size() == 0) {
            for (IndexType nodeId : blk2node[0]) {
                if (blk2node[0].size() == 1) {
                    node_gains[nodeId] += edge_weight;
                } else {
                    node_gains[nodeId] -= edge_weight;
                }
            }
        }

        if (blk2node[1].size() == 1 || blk2node[0].size() == 0) {
            for (IndexType nodeId : blk2node[1]) {
                if (blk2node[1].size() == 1) {
                    node_gains[nodeId] += edge_weight;
                } else {
                    node_gains[nodeId] -= edge_weight;
                }
            }
        }
    }
    return node_gains;
}

template<typename NodeWeightType, typename EdgeWeightType>
EdgeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_gain_of_node(IndexType nodeId) const {
    MoveGainType node_gain = gain_bucket_array.get_value_of(nodeId);
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
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_blk_max_size() const {
    NodeWeightType total_node_weight = hypergraph.get_total_node_weight();
    double half_node_weight = (double) total_node_weight / 2.0;
    return std::floor(half_node_weight * (1 + epsilon));
}

template<typename NodeWeightType, typename EdgeWeightType>
NodeWeightType FMPartitioner<NodeWeightType, EdgeWeightType>::get_blk_min_size() const {
    BlkSizeType total_node_weight = hypergraph.get_total_node_weight();
    double half_node_weight = (double) total_node_weight / 2.0;
    return std::ceil(half_node_weight * (1 - epsilon));
}

template<typename NodeWeightType, typename EdgeWeightType>
std::string FMPartitioner<NodeWeightType, EdgeWeightType>::get_partition_info() {
    std::stringstream ss;
    ss << "cut_size=" << cut_size;
    ss << " blk0=" << blk_sizes[0];
    ss << " blk1=" << blk_sizes[1];
    return ss.str();
}

template<typename NodeWeightType, typename EdgeWeightType>
bool FMPartitioner<NodeWeightType, EdgeWeightType>::is_move_legal(IndexType nodeId, IndexType blkId, std::array<BlkSizeType, 2> tmp_blk_sizes) const {
    
    NodeWeightType node_weight = hypergraph.get_weight_of_node(nodeId);
    IndexType from_blk_id = blkId;
    IndexType to_blk_id = get_opposite_blk_id(from_blk_id);
    tmp_blk_sizes[from_blk_id] -= node_weight;
    tmp_blk_sizes[to_blk_id] += node_weight;

    return tmp_blk_sizes[to_blk_id] <= blk_max_size && tmp_blk_sizes[from_blk_id] >= blk_min_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
void FMPartitioner<NodeWeightType, EdgeWeightType>::move_node(
    IndexType nodeId,
    PartitionResType& tmp_part_res, 
    std::array<BlkSizeType, 2>& tmp_blk_sizes,
    BucketArrayType& tmp_bucket_array
) {
    // move node and update partition states
    NodeWeightType node_weight = hypergraph.get_weight_of_node(nodeId);

    // update node states
    IndexType from_blk_id = tmp_part_res[nodeId];
    tmp_blk_sizes[from_blk_id] -= node_weight;
    IndexType to_blk_id = get_opposite_blk_id(from_blk_id);
    tmp_blk_sizes[to_blk_id] += node_weight;

    // update bucket array
    update_bucket_array(nodeId, tmp_part_res, tmp_bucket_array);

    // update partition results
    tmp_part_res[nodeId] = to_blk_id;
}

template<typename NodeWeightType, typename EdgeWeightType>
void FMPartitioner<NodeWeightType, EdgeWeightType>::update_bucket_array(IndexType nodeId, const PartitionResType& part_res, BucketArrayType& bucket_array) {
    // update the gain bucket array after moving specified node

    IndexType from_blk_id = part_res[nodeId];
    IndexType to_blk_id = get_opposite_blk_id(from_blk_id);

    for (const IndexType& edgeId : hypergraph.get_edges_of_node(nodeId)) {
        EdgeWeightType edge_weight = hypergraph.get_weight_of_edge(edgeId);
        
        std::array<std::vector<IndexType>, 2> blk2node;
        for (IndexType eNodeId : hypergraph.get_nodes_of_edge(edgeId)) {
            if (eNodeId == nodeId) continue;
            IndexType blkId = part_res[eNodeId];
            blk2node[blkId].push_back(eNodeId);
        }

        if (blk2node[to_blk_id].size() == 0) {
            for (IndexType eNodeId : blk2node[from_blk_id]) {
                if (!bucket_array.has_elem(eNodeId)) continue;
                MoveGainType node_gain = bucket_array.get_value_of(eNodeId);
                node_gain += edge_weight;
                bucket_array.update(eNodeId, node_gain);
            }

            if (bucket_array.has_elem(nodeId)) {
                MoveGainType node_gain = bucket_array.get_value_of(nodeId);
                node_gain += 2 * edge_weight;
                bucket_array.update(nodeId, node_gain);                
            }
        }

        if (blk2node[to_blk_id].size() == 1) {
            IndexType eNodeId = blk2node[to_blk_id][0];
            if (bucket_array.has_elem(eNodeId)) {
                MoveGainType node_gain = bucket_array.get_value_of(eNodeId);
                node_gain -= edge_weight;
                bucket_array.update(eNodeId, node_gain);                
            }
        }

        if (blk2node[from_blk_id].size() == 0) {
            for (IndexType eNodeId : blk2node[to_blk_id]) {
                if (!bucket_array.has_elem(eNodeId)) continue;
                MoveGainType node_gain = bucket_array.get_value_of(eNodeId);
                node_gain -= edge_weight;
                bucket_array.update(eNodeId, node_gain);
            }

            if (bucket_array.has_elem(nodeId)) {
                MoveGainType node_gain = bucket_array.get_value_of(nodeId);
                node_gain -= 2 * edge_weight;
                bucket_array.update(nodeId, node_gain);                
            }
        }

        if (blk2node[from_blk_id].size() == 1) {
            IndexType eNodeId = blk2node[from_blk_id][0];
            if (bucket_array.has_elem(eNodeId)) {
                MoveGainType node_gain = bucket_array.get_value_of(eNodeId);
                node_gain += edge_weight;
                bucket_array.update(eNodeId, node_gain);                
            }
        }
    }
}

} // namespace fm_partitioner
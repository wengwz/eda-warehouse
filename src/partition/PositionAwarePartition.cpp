
#include <iostream>
#include <sstream>
#include <random>
#include "PositionAwarePartition.h"

namespace partition {

template<typename NodeWeightType, typename EdgeWeightType>
PositionAwarePartition<NodeWeightType, EdgeWeightType>::PositionAwarePartition(
    int width, int height,
    const std::vector<NodeWeightType>& blk_size_constr,
    const NetlistGraphType& netlist_graph
) {
    // construct partition block graph and netlist graph
    build_blk_graph(width, height, blk_size_constr);
    this->netlist_graph = netlist_graph;

    std::cout << "Netlist Info: " << std::endl << get_netlist_info() << std::endl;
    // init partition
    init_partition();

    std::cout << "Block Info: " << std::endl << get_block_info() << std::endl;
    std::cout << "Netlist Info: " << std::endl << get_netlist_info() << std::endl;
    std::cout << "Initial Partition Info: " << std::endl << get_partition_info() << std::endl;
}

template<typename NodeWeightType, typename EdgeWeightType>
void PositionAwarePartition<NodeWeightType, EdgeWeightType>::run() {

    IndexType iter_count = 0;
    
    std::random_device rd;
    std::default_random_engine rd_gen(rd());

    std::cout << "Start Edge-Based Partition Refinement" << std::endl;
    while (true) {
        std::cout << "Iteration-" << iter_count << ": "<< std::endl;
        std::vector<IndexType> edges;
        for (IndexType edgeId = 0; edgeId < netlist_graph.get_edges_num(); edgeId++) {
            edges.push_back(edgeId);
        }
        std::shuffle(edges.begin(), edges.end(), rd_gen);

        MoveGainType iter_gain = 0;
        for (IndexType edgeId : edges) {
            std::vector<IndexType> edge_blks;
            std::map<IndexType, std::vector<IndexType>> blk2node_map;
            for (IndexType nodeId : netlist_graph.get_nodes_of_edge(edgeId)) {
                IndexType blkId = node2blk_id[nodeId];
                if (blk2node_map.count(blkId) == 0) {
                    edge_blks.push_back(blkId);
                    blk2node_map[blkId] = {nodeId};
                } else {
                    blk2node_map[blkId].push_back(nodeId);
                }
            }
            assert(blk2node_map.size() <= 2);

            if (blk2node_map.size() == 1) continue;
            auto [from_blkId, moves_gain] = get_max_gain_blk(blk2node_map);
            if (from_blkId != 1 && moves_gain > 0) {
                IndexType to_blkId = edge_blks[0] == from_blkId ? edge_blks[1] : edge_blks[0];
                std::cout << "Trial Move nodes of edge-" << edgeId << " from blk-" << from_blkId << " to blk-" << to_blkId << " gain=" << moves_gain << std::endl;
                for (IndexType nodeId : blk2node_map[from_blkId]) {
                    move_node(nodeId, to_blkId);
                }
                iter_gain += moves_gain;
            }
        }

        if (iter_gain == 0) {
            break;
        }
        iter_count++;
    }
    std::cout << "Complete Edge-Based Partition Refinement" << std::endl;

    std::cout << "Start Node-Based Partition Refinement" << std::endl;
    iter_count = 0;
    while (true) {
        std::cout << "Iteration-" << iter_count << ": "<< std::endl;
        MoveGainType iter_gain = 0;
        
        std::vector<IndexType> nodes;
        for (IndexType nodeId = 0; nodeId < netlist_graph.get_nodes_num(); nodeId++) {
            nodes.push_back(nodeId);
        }
        std::shuffle(nodes.begin(), nodes.end(), rd_gen);

        std::vector<IndexType> origin_node2blk_id = node2blk_id;

        for (IndexType nodeId : nodes) {
            auto [max_gain_blk, max_gain] = get_max_gain_blk(nodeId);

            if (max_gain_blk != -1 && max_gain >= 0) {
                std::cout << "Trial Move node-" << nodeId << " from blk-" << node2blk_id[nodeId] << " to blk-" << max_gain_blk << " gain=" << max_gain << std::endl;
                move_node(nodeId, max_gain_blk);
                iter_gain += max_gain;
            }
        }

        if (iter_gain == 0) {
            for (IndexType nodeId = 0; nodeId < netlist_graph.get_nodes_num(); nodeId++) {
                move_node(nodeId, origin_node2blk_id[nodeId]);
            }
        }

        std::cout << "Total Gain of Iteration-" << iter_count << ": " << iter_gain << std::endl;
        std::cout << "Partition Results of Iteration-" << iter_count << ": " << get_partition_info() << std::endl;
        iter_count++;

        if (iter_gain == 0) {
            break;
        }
    }

    std::cout << "Complete Node-Based Partition Refinement" << std::endl;
    //

    check_partition();
}


// private helper functions

template<typename NodeWeightType, typename EdgeWeightType>
void PositionAwarePartition<NodeWeightType, EdgeWeightType>::build_blk_graph(
    int width, int height,
    const std::vector<NodeWeightType>& blk_size_constr
) {

    // partition blocks are assumed to be arranged as array:
    // an example of 3x2 array of partition blocks and its corresponding connectivity:
    // 3 —— 4 —— 5
    // |    |    |
    // 0 —— 1 —— 2

    // add block nodes
    int blk_num = width * height;
    assert(blk_size_constr.size() == blk_num);
    for (int i = 0; i < blk_num; ++i) {
        block_graph.add_node(blk_size_constr[i]);
    }

    // add edge for vertically adjacent block nodes
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < (height - 1); ++y) {
            int adj_blk0 = x + y * width;
            int adj_blk1 = x + (y + 1) * width;
            block_graph.add_edge(1, {adj_blk0, adj_blk1});
        }
    }

    // add edge for horizontally adjacent block nodes
    for (int x = 0; x < (width - 1); ++x) {
        for (int y = 0; y < height; ++y) {
            int adj_blk0 = x + y * width;
            int adj_blk1 = adj_blk0 + 1;
            block_graph.add_edge(1, {adj_blk0, adj_blk1});
        }
    }

    block_graph.build_node2edge_map();

    // build block2neighbor map
    std::cout << "Building block2neighbor map" << std::endl;
    for (IndexType blkId = 0; blkId < block_graph.get_nodes_num(); blkId++) {
        Dist2IndexMap dist2idx_map = get_dist2idx_map(block_graph, blkId, width + height);

        for (int dist = 1; dist < dist2idx_map.size(); ++dist) {
            for (IndexType id : dist2idx_map[dist - 1]) {
                dist2idx_map[dist].push_back(id);
            }
        }
        blk2neighbors_map.push_back(dist2idx_map);
    }
}

template<typename NodeWeightType, typename EdgeWeightType>
void PositionAwarePartition<NodeWeightType, EdgeWeightType>::init_partition() {
    std::cout << "Start Initial Partition" << std::endl;

    CddtNumPriorityBucket cddt_blk_num_pq(netlist_graph.get_nodes_num());
    // init partition results
    cut_size = 0;
    blk_sizes.resize(get_blk_num(), 0);

    node2blk_id.resize(netlist_graph.get_nodes_num(), -1);
    for (IndexType nodeId = 0; nodeId < netlist_graph.get_nodes_num(); nodeId++) {
        node2cddt_blks.push_back(std::vector<bool>(get_blk_num(), true));
    }

    // insert nodes into priority queue
    using pair_type = std::pair<NodeWeightType, IndexType>;
    struct NodeWeightCmp
    {
        bool operator()(const pair_type p1, pair_type p2) {
            return p1.first > p2.first;
        }
    };
    std::priority_queue<pair_type, std::vector<pair_type>, NodeWeightCmp> node_weight_pq;
    for (IndexType nodeId = 0; nodeId < netlist_graph.get_nodes_num(); nodeId++) {
        node_weight_pq.push({netlist_graph.get_weight_of_node(nodeId), nodeId});
    }
    while (!node_weight_pq.empty())
    {
        pair_type node_pair = node_weight_pq.top();
        node_weight_pq.pop();
        cddt_blk_num_pq.insert(node_pair.second, get_blk_num());
    }
    
    // initialize each node to blocks
    IndexType node_count = 0;
    while (!cddt_blk_num_pq.empty()) {
        IndexType nodeId = cddt_blk_num_pq.get_top_id();

        std::cout << std::endl << "Init Partition of node-" << nodeId  << " " << netlist_graph.get_weight_of_node(nodeId) << std::endl;
        std::cout << "Num of Candidates: " << get_cddt_blks_num(node2cddt_blks[nodeId]) << std::endl;
        
        IndexType max_gain_blk = get_initial_blk(nodeId);
        if (max_gain_blk == -1) {
            show_neighbors(nodeId);
        }
        assert(max_gain_blk != -1);

        move_node(nodeId, max_gain_blk); // update partition states

        cddt_blk_num_pq.pop();

        // update cddt_blk_num_pq
        for (IndexType nNodeId : get_neighbors_of_node(nodeId, get_max_dist_of_blk(max_gain_blk))) {
            cddt_blk_num_pq.update(nNodeId, get_cddt_blks_num(node2cddt_blks[nNodeId]));
        }
        
        std::cout << "Node-" << nodeId << " is initialized to block-" << max_gain_blk << std::endl;
        std::cout << "Partition Info: " << get_partition_info() << std::endl;
        std::cout << "Num of Nodes Initialized: " << ++node_count << std::endl << std::endl;
    }

    std::cout << "Complete Initial Partition" << std::endl;
}

template<typename NodeWeightType, typename EdgeWeightType>
typename PositionAwarePartition<NodeWeightType, EdgeWeightType>::BlockGainPair 
PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_max_gain_blk(IndexType nodeId) {
    //std::cout << "Start finding max gain block for node-" << nodeId << std::endl;
    IndexType max_gain = 0;
    IndexType max_gain_blk = -1;

    // calculate original num of cut related with nodeId
    IndexType orig_blk_id = node2blk_id[nodeId];
    CutSizeType orig_cut_size = get_cut_size_of_node(nodeId);

    const std::vector<bool>& cddt_blks = node2cddt_blks[nodeId];
    for (IndexType blkId = 0; blkId < get_blk_num(); ++blkId) {
        if (!cddt_blks[blkId] || blkId == get_blk_of_node(nodeId)) {
            continue;
        }

        NodeWeightType new_blk_size = blk_sizes[blkId] + netlist_graph.get_weight_of_node(nodeId);
        NodeWeightType blk_size_limit = block_graph.get_weight_of_node(blkId);

        if (new_blk_size > blk_size_limit) {
            continue; // avoid overflow of block
        }

        // trail move nodeId to blkId and get new cut size
        node2blk_id[nodeId] = blkId;
        CutSizeType trail_cut_size = get_cut_size_of_node(nodeId);
        // recover partition result
        node2blk_id[nodeId] = orig_blk_id;
        CutSizeType move_gain = orig_cut_size - trail_cut_size;

        if (max_gain_blk == -1 || move_gain > max_gain) {
            max_gain = move_gain;
            max_gain_blk = blkId;
        }
    }
    return {max_gain_blk, max_gain};
}

template<typename NodeWeightType, typename EdgeWeightType>
typename PositionAwarePartition<NodeWeightType, EdgeWeightType>::BlockGainPair
PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_max_gain_blk(std::map<IndexType, std::vector<IndexType>>& blk2node_map) {
    // get block to be moved with maximum gain for a cut edge
    IndexType max_gain_blk = -1;
    CutSizeType max_gain = 0;

    std::vector<IndexType> edge_blks;
    for (const auto& [blkId, nodes]: blk2node_map) {
        edge_blks.push_back(blkId);
    }
    assert(edge_blks.size() == 2);

    for (int i = 0; i < edge_blks.size(); i++) {
        IndexType from_blkId = edge_blks[i];
        IndexType to_blkId = edge_blks[1 - i];

        const std::vector<IndexType>& from_blk_nodes = blk2node_map[from_blkId];

        // check if move from_blkId to to_blkId is legal
        bool all_has_cddt_blk = true;
        BlkSizeType new_blk_size = blk_sizes[to_blkId];
        for (IndexType nodeId : from_blk_nodes) {
            if (!node2cddt_blks[nodeId][to_blkId]) {
                all_has_cddt_blk = false;
            }
            new_blk_size += netlist_graph.get_weight_of_node(nodeId);
        }
        BlkSizeType blk_size_limit = block_graph.get_weight_of_node(to_blkId);
        if (new_blk_size > blk_size_limit || !all_has_cddt_blk) {
            continue; // avoid overflow of block
        }

        // calculate gain of moving from_blkId to to_blkId
        MoveGainType moves_gain = 0;
        for (IndexType nodeId : from_blk_nodes) {
            CutSizeType orig_cut_size = get_cut_size_of_node(nodeId);
            node2blk_id[nodeId] = to_blkId;
            CutSizeType new_cut_size = get_cut_size_of_node(nodeId);
            moves_gain += orig_cut_size - new_cut_size;
        }

        // recover partition result
        for (IndexType nodeId : from_blk_nodes) {
            node2blk_id[nodeId] = from_blkId;
        }

        if (max_gain_blk == -1 || moves_gain > max_gain) {
            max_gain_blk = from_blkId;
            max_gain = moves_gain;
        }
    }

    return {max_gain_blk, max_gain};
}

template<typename NodeWeightType, typename EdgeWeightType>
IndexType PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_initial_blk(IndexType nodeId) {
    // get initial block for node
    // return the block with minimum size
    NodeWeightType min_blk_size;
    IndexType min_blk_id = -1;

    const std::vector<bool>& cddt_blks = node2cddt_blks[nodeId];
    for (IndexType blkId = 0; blkId < get_blk_num(); ++blkId) {
        if (cddt_blks[blkId] == false || blkId == get_blk_of_node(nodeId)) {
            continue;
        }

        // check if move nodeId to blkId results in block overflow
        BlkSizeType blk_size_limit = block_graph.get_weight_of_node(blkId);
        BlkSizeType new_blk_size = blk_sizes[blkId] + netlist_graph.get_weight_of_node(nodeId);
        if (new_blk_size > blk_size_limit) {
            continue; // avoid overflow of block
        }

        // check if this move results in no candidate blks for neighboring nodes
        Dist2IndexMap dist2neighbors = get_dist2idx_map(netlist_graph, nodeId, get_max_dist_of_blk(blkId));
        bool all_has_cddt_blk = true;
        for (int dist = 1; dist < dist2neighbors.size(); dist++) {

            for (IndexType nNodeId : dist2neighbors[dist]) {
                std::vector<bool> new_cddt_blks = merge_cddt_blks(get_cddt_blks(blkId, dist), node2cddt_blks[nNodeId]);

                if (!has_cddt_blks(new_cddt_blks)) {
                    all_has_cddt_blk = false;
                    std::cout << "No candidate blks for node-" << nNodeId << "if move to blk-" << blkId << std::endl;
                    break;
                }
            }

            if (!all_has_cddt_blk) break;
        }

        if (all_has_cddt_blk) {
            if (min_blk_id == -1 || new_blk_size < min_blk_size) {
                min_blk_size = new_blk_size;
                min_blk_id = blkId;
            }
        }
    }

    return min_blk_id;
}

template<typename NodeWeightType, typename EdgeWeightType>
typename PositionAwarePartition<NodeWeightType, EdgeWeightType>::CutSizeType 
PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_cut_size_of_node(IndexType nodeId) const {
    return get_cut_size_of_node(nodeId, node2blk_id);
}

template<typename NodeWeightType, typename EdgeWeightType>
typename PositionAwarePartition<NodeWeightType, EdgeWeightType>::CutSizeType
PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_cut_size_of_node(IndexType nodeId, const std::vector<IndexType>& trial_node2blk) const {
    CutSizeType node_cut_size = 0;
    for (IndexType edgeId : netlist_graph.get_edges_of_node(nodeId)) {
        std::set<IndexType> edge2blks;

        for (IndexType eNodeId : netlist_graph.get_nodes_of_edge(edgeId)) {
            if (trial_node2blk[eNodeId] == -1) continue;
            edge2blks.insert(trial_node2blk[eNodeId]);
        }

        if (edge2blks.size() > 1) {
            node_cut_size += netlist_graph.get_weight_of_edge(edgeId);
        }
    }

    return node_cut_size;
}

template<typename NodeWeightType, typename EdgeWeightType>
bool PositionAwarePartition<NodeWeightType, EdgeWeightType>::move_node(IndexType nodeId, IndexType new_blk) {
    // move nodeId to new_blk and update partition states
    if (node2blk_id[nodeId] == new_blk) {
        return true;
    }

    IndexType orig_blk = node2blk_id[nodeId];
    CutSizeType orig_cut_size = get_cut_size_of_node(nodeId);

    //
    if (!node2cddt_blks[nodeId][new_blk]) {
        return false;
    }

    node2blk_id[nodeId] = new_blk;

    // update candidate block of neighbors
    int max_dist = get_max_dist_of_blk(new_blk);
    Dist2IndexMap dist2node_map = get_dist2idx_map(netlist_graph, nodeId, max_dist);
    
    for (int dist = 1; dist < dist2node_map.size(); dist++) {
        std::vector<bool> new_cddt_blks(get_blk_num(), false);
        for (IndexType blk : blk2neighbors_map[new_blk][dist]) {
            new_cddt_blks[blk] = true;
        }

        for (IndexType nNodeId : dist2node_map[dist]) {
            node2cddt_blks[nNodeId] = merge_cddt_blks(node2cddt_blks[nNodeId], new_cddt_blks);
        }
    }

    // update cut size
    CutSizeType new_cut_size = get_cut_size_of_node(nodeId);
    cut_size += new_cut_size - orig_cut_size;

    // update block size
    NodeWeightType node_weight = netlist_graph.get_weight_of_node(nodeId);
    blk_sizes[new_blk] += node_weight;
    if (orig_blk >= 0) {
        blk_sizes[orig_blk] -= node_weight;
    }
    
    return true;
}

template<typename NodeWeightType, typename EdgeWeightType>
const std::vector<IndexType>& PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_neighbors_of_blk(IndexType blkId, int dist) const {
    if (dist >= blk2neighbors_map[blkId].size()) {
        return blk2neighbors_map[blkId].back();
    }
    return blk2neighbors_map[blkId][dist];
}

template<typename NodeWeightType, typename EdgeWeightType>
std::vector<IndexType> PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_neighbors_of_node(IndexType nodeId, int dist) {
    Dist2IndexMap dist2node_map = get_dist2idx_map(netlist_graph, nodeId, dist);
    std::vector<IndexType> neighbor_nodes;

    for (int i = 1; i < dist2node_map.size(); i++) {
        for (IndexType id : dist2node_map[i]) {
            neighbor_nodes.push_back(id);
        }
    }

    return neighbor_nodes;
}

template<typename NodeWeightType, typename EdgeWeightType>
typename PositionAwarePartition<NodeWeightType, EdgeWeightType>::Dist2IndexMap
PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_dist2idx_map(const NetlistGraphType& graph, IndexType nodeId, int maxDist) {
    Dist2IndexMap dist2node_map;

    std::map<IndexType, int> node2dist_map;

    std::queue<IndexType> node_queue;
    node_queue.push(nodeId);
    node2dist_map[nodeId] = 0;
    dist2node_map.push_back({nodeId});

    while (!node_queue.empty()) {
        IndexType cur_node = node_queue.front();
        node_queue.pop();

        for (IndexType id : graph.get_neighbors_of_node(cur_node)) {
            if (node2dist_map.count(id) == 0) {
                int dist = node2dist_map[cur_node] + 1;
                if (dist > maxDist) continue;

                node_queue.push(id);
                node2dist_map[id] = dist;
                
                dist2node_map.resize(dist + 1, std::vector<IndexType>());
                dist2node_map[dist].push_back(id);
            }
        }
    }
    return dist2node_map;
}


template<typename NodeWeightType, typename EdgeWeightType>
std::vector<bool> PositionAwarePartition<NodeWeightType, EdgeWeightType>::merge_cddt_blks(const std::vector<bool>& blks1, const std::vector<bool>& blks2) {
    assert(blks1.size() == blks2.size());
    std::vector<bool> result;
    for (int i = 0; i < blks1.size(); i++) {
        result.push_back(blks1[i] && blks2[i]);
    }
    return result;
}

template<typename NodeWeightType, typename EdgeWeightType>
bool PositionAwarePartition<NodeWeightType, EdgeWeightType>::has_cddt_blks(const std::vector<bool>& cddt_blks) const {
    return std::any_of(cddt_blks.begin(), cddt_blks.end(), [](bool blk) {return blk;});
}

template<typename NodeWeightType, typename EdgeWeightType>
IndexType PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_cddt_blks_num(const std::vector<bool>& cddt_blks) const {
    return std::count(cddt_blks.begin(), cddt_blks.end(), true);
}

template<typename NodeWeightType, typename EdgeWeightType>
std::vector<bool> PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_cddt_blks(IndexType blkId, int dist) const {
    std::vector<bool> cddt_blks(get_blk_num(), false);
    for (IndexType blkId : blk2neighbors_map[blkId][dist]) {
        cddt_blks[blkId] = true;
    }
    return cddt_blks;
}

template<typename NodeWeightType, typename EdgeWeightType>
void PositionAwarePartition<NodeWeightType, EdgeWeightType>::show_neighbors(IndexType nodeId) {

    Dist2IndexMap dist2idx_map = get_dist2idx_map(netlist_graph, nodeId, 3);
    for (int i = 0; i < dist2idx_map.size(); i++) {
        std::cout << "Dist-" << i << ": ";
        for (IndexType id : dist2idx_map[i]) {
            std::cout << id << ":" << node2blk_id[id] << " ";
        }
        std::cout << std::endl;
    }
}


template<typename NodeWeightType, typename EdgeWeightType>
void PositionAwarePartition<NodeWeightType, EdgeWeightType>::check_partition() {
    // check the legality of partition results
    std::cout << "Check Partition Result:" << std::endl;
    CutSizeType check_cut_size = 0;
    std::map<IndexType, IndexType> edge_fanout2num_map;
    for (IndexType edgeId = 0; edgeId < netlist_graph.get_edges_num(); edgeId++) {
        std::set<IndexType> edge2blks;
        for (IndexType nodeId : netlist_graph.get_nodes_of_edge(edgeId)) {
            assert(node2blk_id[nodeId] != -1);
            edge2blks.insert(node2blk_id[nodeId]);
        }

        if (edge2blks.size() > 1) {
            check_cut_size += netlist_graph.get_weight_of_edge(edgeId);
        }

        if (edge_fanout2num_map.count(edge2blks.size()) == 0) {
            edge_fanout2num_map[edge2blks.size()] = 1;
        } else {
            edge_fanout2num_map[edge2blks.size()]++;
        }
    }
    std::cout << "Check Cut Size: " << check_cut_size << std::endl;

    for (auto& pair : edge_fanout2num_map) {
        std::cout << "fanout-" << pair.first << ": " << pair.second << std::endl;
    }

    // check block sizes
    std::vector<BlkSizeType> check_blk_sizes(get_blk_num(), 0);
    for (IndexType nodeId = 0; nodeId < netlist_graph.get_nodes_num(); nodeId++) {
        IndexType blkId = node2blk_id[nodeId];
        assert(blkId != -1);
        check_blk_sizes[blkId] += netlist_graph.get_weight_of_node(nodeId);
    }

    std::cout << "Check Block Sizes: ";
    for (IndexType i = 0; i < get_blk_num(); i++) {
        std::cout << "blk-" << i << " " << check_blk_sizes[i] << " ";
    }
    std::cout << std::endl;

    std::map<IndexType, IndexType> cddt_num2num_map;
    for (IndexType nodeId = 0; nodeId < netlist_graph.get_nodes_num(); nodeId++) {
        IndexType cddt_num = get_cddt_blks_num(node2cddt_blks[nodeId]);
        if (cddt_num2num_map.count(cddt_num) == 0) {
            cddt_num2num_map[cddt_num] = 1;
        } else {
            cddt_num2num_map[cddt_num]++;
        }
    }

    std::cout << "Check Num of Candidate Blocks: ";
    for (auto& pair : cddt_num2num_map) {
        std::cout << "cddt-" << pair.first << ": " << pair.second << " ";
    }
    std::cout << std::endl;
}


template<typename NodeWeightType, typename EdgeWeightType>
std::string PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_partition_info() const {
    std::stringstream info;
    info << "Cut Size: " << cut_size << std::endl;
    info << "Block Sizes: ";
    for (IndexType i = 0; i < get_blk_num(); i++) {
        info << "Blk-" << i << "=" << blk_sizes[i] << " ";
    }
    return info.str();
}

template<typename NodeWeightType, typename EdgeWeightType>
std::string PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_block_info() const {
    std::stringstream info;

    info << "Num of Blocks: " << block_graph.get_nodes_num() << std::endl;
    info << "Block Size Limitation: ";
    for (IndexType i = 0; i < block_graph.get_nodes_num(); i++) {
        info << "Blk-" << i << "=" << block_graph.get_weight_of_node(i) << " ";
    }
    info << std::endl << "Neighbors of Blocks: " << std::endl;
    for (IndexType blkId = 0; blkId < block_graph.get_nodes_num(); blkId++) {
        info << "block-" << blkId << ":" << std::endl;
        for (IndexType dist = 0; dist <= get_max_dist_of_blk(blkId); dist++) {
            info << "  dist <= " << dist << ": ";
            for (IndexType neighborId : get_neighbors_of_blk(blkId, dist)) {
                info << neighborId << " ";
            }
            info << std::endl;
        }
    }
    return info.str();
}

template<typename NodeWeightType, typename EdgeWeightType>
std::string PositionAwarePartition<NodeWeightType, EdgeWeightType>::get_netlist_info() const {
    std::stringstream info;

    info << "Num of Nodes: " << netlist_graph.get_nodes_num() << " " << "Total Weight: " << netlist_graph.get_total_node_weight() << " " << std::endl;
    info << "Num of Edges: " << netlist_graph.get_edges_num() << " " << "Total Weight: " << netlist_graph.get_total_edge_weight() << " ";
    return info.str();
}

} // namespace partition
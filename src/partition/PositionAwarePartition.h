
#ifndef POISITION_AWARE_PARTITION_H
#define POISITION_AWARE_PARTITION_H

#include <vector>
#include <queue>
#include <map>
#include <set>

#include "hypergraph.h"
#include "bucket_array.h"

using namespace hypergraph;
namespace partition {

template<typename NodeWeightType, typename EdgeWeightType>
class PositionAwarePartition {
    using NetlistGraphType = HyperGraph<NodeWeightType, EdgeWeightType>;
    using BlockGraphType = NetlistGraphType;

    using Dist2IndexMap = std::vector<std::vector<IndexType>>;
    using Blk2NeighborsMap = std::vector<Dist2IndexMap>;

    using MoveGainType = EdgeWeightType;
    using CutSizeType = EdgeWeightType;
    using BlkSizeType = NodeWeightType;

    using CddtNumPriorityBucket = MinPriorityBucketArray<IndexType, IndexType>;

public:
    PositionAwarePartition(
        int width, int height,
        const std::vector<NodeWeightType>& blk_size_constr,
        const NetlistGraphType& netlist_graph
    );

    void run(); // kernel function of partitioning


public:
    // getters
    IndexType get_blk_num() const {
        return block_graph.get_nodes_num();
    }

    IndexType get_node_num() const {
        return netlist_graph.get_nodes_num();
    }
    
    IndexType get_blk_of_node(IndexType nodeId) const {
        return node2blk_id[nodeId];
    }

private:
    void build_blk_graph(int width, int height, const std::vector<NodeWeightType>& blk_size_constr);
    void init_partition();

    IndexType get_max_gain_blk(IndexType nodeId, const std::vector<bool>& cddt_blks);
    CutSizeType get_cut_size_of_node(IndexType nodeId) const;
    CutSizeType get_cut_size_of_node(IndexType nodeId, const std::vector<IndexType>& node2blk) const;
    //std::vector<IndexType> move_node(IndexType nodeId, IndexType new_blk);
    bool move_node(IndexType nodeId, IndexType new_blk);

    const std::vector<IndexType>& get_neighbors_of_blk(IndexType blkId, int dist) const;
    std::vector<IndexType> get_neighbors_of_node(IndexType nodeId, int dist);
    int get_max_dist_of_blk(IndexType blkId) const {
        return blk2neighbors_map[blkId].size() - 1;
    }

    Dist2IndexMap get_dist2idx_map(const NetlistGraphType& graph, IndexType id, int maxDist);

    std::vector<bool> merge_cddt_blks(const std::vector<bool>& blks1, const std::vector<bool>& blks2);
    bool has_cddt_blks(const std::vector<bool>& cddt_blks) const;
    IndexType get_cddt_blks_num(const std::vector<bool>& cddt_blks) const;

    std::vector<bool> get_cddt_blks(IndexType blkId, int dist) const;
    std::vector<bool> get_cddt_blks(IndexType nodeId) const;

    // debug functions
    void check_partition();
    void show_neighbors(IndexType nodeId);
    std::string get_partition_info() const;
    std::string get_block_info() const;
    std::string get_netlist_info() const;

private:

    // hypergraph of netlist and blocks
    NetlistGraphType netlist_graph;
    BlockGraphType block_graph;
    Blk2NeighborsMap blk2neighbors_map;

    // partition results
    //// edge states
    CutSizeType cut_size;

    //// block states
    std::vector<BlkSizeType> blk_sizes;
    
    //// node states
    std::vector<IndexType> node2blk_id;
    std::vector<std::vector<bool>> node2cddt_blks;
};

template class PositionAwarePartition<int, int>;

}

#endif
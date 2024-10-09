
#ifndef FM_PARTITIONER_H
#define FM_PARTITIONER_H

#include "hypergraph.h"
#include "bucket_array.h"

namespace fm_partitioner {
using namespace hypergraph;

// 2-Way Fiduccia-Mattheyses Partitioner
template<typename NodeWeightType, typename EdgeWeightType>
class FMPartitioner {
    using HyperGraphType = HyperGraph<NodeWeightType, EdgeWeightType>;
    using MoveGainType = EdgeWeightType;
    using BlkSizeType = NodeWeightType;
    using BucketArrayType = GainBucketArray<MoveGainType, IndexType>;
    using PartitionResType = std::vector<IndexType>;

public:
    FMPartitioner(const HyperGraphType& hypergraph, double epsilon) 
    : hypergraph(hypergraph), gain_bucket_array(hypergraph.get_nodes_num()), epsilon(epsilon) {
        init_partition();
    };

    FMPartitioner(
        const std::vector<NodeWeightType>& node_weights,
        const std::vector<EdgeWeightType>& edge_weights,
        const std::vector<std::vector<IndexType>>& edge2node, 
        double epsilon
    );

    void operator()() { run() };

    void run(); // kernel function of partitioning

public:
    EdgeWeightType get_cutsize() const;
    BlkSizeType get_size_of_blk(IndexType blkId) const;
    std::vector<IndexType> get_nodes_of_blk(IndexType blkId) const;
    IndexType get_blk_id_of_node(IndexType node) const;
    const std::vector<IndexType>& get_partition_result() const;

private:
    void init_partition();

    EdgeWeightType get_cutsize(const std::vector<IndexType> node_blk_ids) const;
    
    MoveGainType get_gain_of_node(IndexType nodeId) const;
    std::vector<MoveGainType> get_gain_of_nodes(const std::vector<IndexType>& tmp_node_blk_ids) const;
    
    BlkSizeType get_size_of_blk(IndexType blkId, const std::vector<IndexType>& tmp_node_blk_ids) const;
    BlkSizeType get_blk_max_size() const;
    BlkSizeType get_blk_min_size() const;
    
    inline IndexType get_opposite_blk_id(IndexType blkId) const { return (blkId == 0) ? 1 : 0;}

    bool is_move_legal(IndexType nodeId, IndexType blkId, std::array<BlkSizeType, 2> tmp_blk_sizes) const;

    std::string get_partition_info();

    void move_node(
        IndexType nodeId,
        PartitionResType& tmp_part_res, 
        std::array<BlkSizeType, 2>& tmp_blk_sizes,
        BucketArrayType& tmp_gain_bucket_array
    );

    void update_bucket_array(IndexType nodeId, const PartitionResType& part_res, BucketArrayType& buckt_array);

private:
    // partition constraints
    double epsilon;
    BlkSizeType blk_max_size; // The maximum size of a block
    BlkSizeType blk_min_size; // The minimum size of a block

    // target hypergraph
    HyperGraph<NodeWeightType, EdgeWeightType> hypergraph;
    
    // partition results
    //// node states
    std::vector<IndexType> node_blk_ids;
    std::array<BlkSizeType, 2> blk_sizes;
    //// edge states
    EdgeWeightType cut_size;
    //// gain bucket array
    BucketArrayType gain_bucket_array;
};

template class FMPartitioner<int, int>;
}

#endif
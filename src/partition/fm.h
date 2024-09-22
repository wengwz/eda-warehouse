
#ifndef FM_PARTITIONER_H
#define FM_PARTITIONER_H

#include "hypergraph.h"

namespace fm_partitioner {
using namespace hypergraph;

// 2-Way Fiduccia-Mattheyses Partitioner
template<typename NodeWeightType, typename EdgeWeightType>
class FMPartitioner {
    using HyperGraphType = HyperGraph<NodeWeightType, EdgeWeightType>;
    using MoveGainType = EdgeWeightType;
    using BlkSizeType = NodeWeightType;

public:
    FMPartitioner(const HyperGraphType& hypergraph, double epsilon) 
    : hypergraph(hypergraph), epsilon(epsilon) {};

    FMPartitioner(
        const std::vector<NodeWeightType>& node_weights,
        const std::vector<EdgeWeightType>& edge_weights,
        const std::vector<std::vector<IndexType>>& edge2node, 
        double epsilon
    );

    void operator()(); // main function of partitioning

public:
    EdgeWeightType get_cutsize() const;
    NodeWeightType get_size_of_blk(IndexType blkId) const;
    std::vector<IndexType> get_nodes_of_blk(IndexType blkId) const;
    IndexType get_blk_id_of_node(IndexType node) const;
    const std::vector<IndexType>& get_blk_ids() const;

private:
    void init_partition();
    MoveGainType get_gain_of_node(IndexType nodeId, const std::vector<IndexType>& tmp_node_blk_ids) const;
    
    BlkSizeType get_size_of_blk(IndexType blkId, const std::vector<IndexType>& tmp_node_blk_ids) const;
    BlkSizeType get_max_blk_size() const;
    BlkSizeType get_min_blk_size() const;
    
    inline IndexType get_opposite_blk_id(IndexType blkId) const { return (blkId == 0) ? 1 : 0;}
    std::string get_partition_info(EdgeWeightType cut_size, const std::array<NodeWeightType, 2>& blk_sizes);

private:
    double epsilon;
    HyperGraph<NodeWeightType, EdgeWeightType> hypergraph;
    std::vector<IndexType> node_blk_ids;

};

template class FMPartitioner<int, int>;
}

#endif
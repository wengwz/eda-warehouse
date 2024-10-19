
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "PositionAwarePartition.h"

struct ResourceTypeUtil
{
    int FF = 0;
    int CARRY = 0;
    int MISCS = 0;
    int LUT = 0;
};

void from_json(const json& j, ResourceTypeUtil& resourceTypeUtil) {
    if (j.contains("FF")) resourceTypeUtil.FF = j["FF"];
    if (j.contains("CARRY")) resourceTypeUtil.CARRY = j["CARRY"];
    if (j.contains("MISCS")) resourceTypeUtil.MISCS = j["MISCS"];
    if (j.contains("LUT")) resourceTypeUtil.LUT = j["LUT"];
}

struct PartitionGroup
{
    int id;
    int primCellNum;
    ResourceTypeUtil resourceTypeUtil;
    std::vector<std::string> groupCellNames;
};

void from_json(const json& j, PartitionGroup& partitionGroup) {
    j.at("id").get_to(partitionGroup.id);
    j.at("primCellNum").get_to(partitionGroup.primCellNum);
    j.at("resourceTypeUtil").get_to(partitionGroup.resourceTypeUtil);
    j.at("groupCellNames").get_to(partitionGroup.groupCellNames);
}

struct PartitionEdge
{
    int id;
    int primCellNum;
    int weight;
    int degree;
    std::vector<int> incidentGroupIds;
};

void from_json(const json& j, PartitionEdge& partitionEdge) {
    partitionEdge.id = j["id"];
    partitionEdge.primCellNum = j["primCellNum"];
    partitionEdge.weight = j["weight"];
    partitionEdge.degree = j["degree"];
    j.at("incidentGroupIds").get_to(partitionEdge.incidentGroupIds);
}

int main(int argc, char** argv) {
    assert(argc == 3);
    std::string input_json_path = argv[1];
    std::string output_json_path = argv[2];

    std::ifstream input_json_file(input_json_path);
    assert(input_json_file.is_open());

    json input_json;
    input_json_file >> input_json;

    // parse circuit netlist info
    int nodes_num = input_json["totalGroupNum"];
    int edges_num = input_json["totalEdgeNum"];
    HyperGraph<int, int> netlist_graph;
    //// parse node information
    std::cout << "Start parsing node information" << std::endl;
    std::vector<PartitionGroup> partitionGroups;
    input_json.at("partitionGroups").get_to(partitionGroups);
    
    for (int i = 0; i < partitionGroups.size(); ++i) {
        assert(i == partitionGroups[i].id);
        // use the amount of LUTs as weight of nodes
        netlist_graph.add_node(partitionGroups[i].resourceTypeUtil.LUT);
    }

    std::cout << "Start parsing edge information" << std::endl;
    //// parse edge information
    std::vector<PartitionEdge> partitionEdges;
    input_json.at("partitionEdges").get_to(partitionEdges);
    for (int i = 0; i < partitionEdges.size(); i++) {
        netlist_graph.add_edge(partitionEdges[i].weight, partitionEdges[i].incidentGroupIds);
    }
    netlist_graph.build_node2edge_map();


    // parse block size constraints
    int grid_width = input_json["gridWidth"];
    int grid_height = input_json["gridHeight"];
    std::vector<int> grid_limits;
    input_json.at("gridLimits").get_to(grid_limits);
    assert(grid_limits.size() == grid_width * grid_height);
    
    
    // run partitioning
    partition::PositionAwarePartition<int, int> partition(
        grid_width, grid_height, grid_limits, netlist_graph
    );
    partition.run();
}
#ifndef LBNODES_H
#define LBNODES_H

#include <vector>
#include <numeric>
#include "LBglobal.h"
#include "LBlatticetypes.h"

/*
 * Nodes: -contains additional informatino about a node, eg. rank, node type etc.
 *        -class Grid will still contain the cartesion indices of the nodes and
 *         the list of node-numbers of the nodes in a node's neighborhood.
 *
 * Work in progress: moving node-information from class Grid to class Nodes
*/

class Nodes
{
public:
    Nodes(const int nNodes): nNodes_(nNodes), rank_(nNodes) {}
private:
    int nNodes_;
    std::vector<int> rank_;
};

#endif // LBNODES_H

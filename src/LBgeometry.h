#ifndef LBGEOMETRY_H
#define LBGEOMETRY_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBgrid.h"
#include "LBnodes.h"
#include "LBboundary.h"
#include "LBhalfwaybb.h"

template<typename DXQY>
std::vector<int> findBulkNodes(const Nodes<DXQY> &nodes)
// makeBulkNodes : make a list of bulk node labels. Here we assume that all
//  fluid nodes are also bulk nodes.
{
    std::vector<int> bulkNodes;
    for (int n = 1; n < nodes.size(); ++n)
        if (nodes.isFluid(n) && nodes.isMyRank(n))
            bulkNodes.push_back(n);
    return bulkNodes;
}


template<typename DXQY>
std::vector<int> findSolidBndNodes(const Nodes<DXQY> &nodes)
{
    std::vector<int> ret; // List of node numbers of all solid boundary nodes for myRank process
    for (int n = 1; n < nodes.size(); n++) { // Loop over all grid nodes excpet the default node (node number = 0)
        if (nodes.isSolidBoundary(n))  ret.push_back(n);
    }
    return ret;
}


template<typename DXQY>
const std::vector<int> findFluidBndNodes(const Nodes<DXQY> &nodes)
{
    std::vector<int> ret; // List of node numbers to all fluid boundary nodes for myRank process
    for (int n = 1; n < nodes.size(); n++) { // Loop over all grid nodes excpet the default node (node number = 0)
        if (nodes.isFluidBoundary(n))  ret.push_back(n);
    }
    return ret;
}


// template <typename DXQY>
template <template <class> class T,  typename DXQY>
T<DXQY> makeFluidBoundary(const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
/* Sets up Bounce back bounray Boundary object by classifying and adding boundary nodes
 *
 * myRank    : rank of the current process
 * grid      : object of the Grid class
 */

{
    std::vector<int> bndNodes = findFluidBndNodes(nodes);  // Find list of boundary nodes
    T<DXQY> bnd(bndNodes.size());  // Set size of the boundary node object

    for (auto nodeNo: bndNodes) {  // For all boundary nodes
        int nBeta = 0, beta[DXQY::nDirPairs_];
        int nGamma = 0, gamma[DXQY::nDirPairs_];
        int nDelta = 0, delta[DXQY::nDirPairs_];

        for (int q = 0; q < DXQY::nDirPairs_; ++q) {
            int neig_q = grid.neighbor(q, nodeNo);
            int neig_q_rev = grid.neighbor(bnd.dirRev(q), nodeNo);

            if (nodes.isFluid(neig_q)) { // FLUID
                if (nodes.isFluid(neig_q_rev)) { // FLUID :GAMMA
                    gamma[nGamma] = q;
                    nGamma += 1;
                }
                else { // SOLID :BETA (q) and BETA_HAT (qRev)
                    beta[nBeta] = q;
                    nBeta += 1;
                }
            }
            else { // SOLID
                if (nodes.isFluid(neig_q_rev)) { // FLUID :BETA (qDir) and BETA_HAT (q)
                    beta[nBeta] = bnd.dirRev(q);
                    nBeta += 1;
                }
                else { // SOLID: DELTA
                    delta[nDelta] = q;
                    nDelta += 1;
                }
            }
        } // END FOR DIR PAIRS
        bnd.addNode(nodeNo, nBeta, beta, nGamma, gamma, nDelta, delta);
    } // For all boundary nodes
    return bnd;
}


#endif // LBGEOMETRY_H

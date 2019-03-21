#ifndef LBBNDMPI_H
#define LBBNDMPI_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"
#include "Input.h"

/********************************************************* PROTOCOL
 * We assume the existence of three types of input-files:
 * (Find better file extension)
 *
 *  1) rank.mpi : Geometry with process partitioning and solids
 *                Solids : 0
 *                Processor rank: 1 + rank (so that rank 0 is 1 etc.)
 *
 *  2) labels_global.mpi : Geometry with local node numbers
 *
 *  3) labels_rank#.mpi : Node number for process # including
 *                        ghost node and boundary-node labels
 *
 *
 * ************* INIT MPI RUN ****************************
 *  ---- How to setup the send and recive protocols -----
 *
 *  General: setup a list that links the current mpi-boundaries
 *           node to the node that it represetns in the
 *           neighboring process.
 *
 *  Example. Scalar from rank a to rank b
 *
 *   node position notation : r = [x0, x1, ...]:
 *
 *     1) Go through all nodes in labels_ranka
 *        and check which are connected to rank b
 *        That is:  rank.mpi[r] == b.
 *        Then we have the number of nodes in the mpi-boundary
 *        linked to b.
 *        Send this number to proc b. Now be knows how much data to
 *        receive.
 *
 *      2) List protocols: (Scalars)
 *        'Rank a node-list '[labels_ranka[r], ... ]  (Recive list)
 *        'Rank b node-list' [labels_global[r], ... ]  (Send list)
 *         Send 'Rank b node-list' to rank b proc. and this gives
 *         the order for the data communcation.
 *         We know now that the first element in the sent array
 *         shold now be assign to the node-label given by
 *         the first entry in 'Rank a node-list'
 *
 *
 *      3)  Fill list communicated by mpi list.
 *          Rank b: Has 'Rank b node-list' []
 *          Send list scalar-mpi-b = [rho('Rank b node-list'[0]), rho('Rank b node-list'[1]), ....]
 *          to rank a rho(Rank a node-list[0]) = scalar-mpi-b[0], etc.
 *
 *********************************************************/


// READ FILES
template <typename DXQY>
void setupBndMpi(MpiFile<DXQY> &localFile, MpiFile<DXQY> &globalFile, MpiFile<DXQY> &rankFile, const int rank, Grid<DXQY> &grid)
{
    std::vector<int> rankNo; // Rank of adjacent processors
    std::vector<int> nNodes; // Number of nodes in the boundary to each adjacent processor.
    std::vector<bool> nodeAdded(grid.num_nodes(), false);
    for (int pos=0; pos < static_cast<int>(localFile.size()); ++pos) {
        int localLabel, globalLabel;
        int nodeRank;

        rankFile.getVal(nodeRank);
        localFile.getVal(localLabel);
        globalFile.getVal(globalLabel);

        if ( (nodeRank != 0) && (localLabel != 0) ) { // Not a solid and not a local default node
            nodeRank -= 1; // Get the actual rank number
            if (rank != nodeRank) {
                std::vector<int>::iterator iter = std::find(rankNo.begin(), rankNo.end(), nodeRank);
                // If the rank is not registered already: add new rankNo and nNodes element
                if ( iter == rankNo.end() ) {
                   rankNo.push_back(nodeRank);
                   nNodes.push_back(1);
                } else if (!nodeAdded[static_cast<std::size_t>(localLabel)])
                { // If it is registered: increase nNodes for '
                    std::size_t index = static_cast<std::size_t>(std::distance(rankNo.begin(), iter));
                    nNodes[index] += 1;
                    nodeAdded[static_cast<std::size_t>(localLabel)] = true;
                    // std::cout << "local label = " << localLabel << std::endl;
                }
            } // End if it is an mpi boundary site
        } // End if site is Not a solid and not a local default node
    } // End for-loop reading map values
    for (std::size_t n = 0; n < rankNo.size(); ++n) {
        std::cout << "rank = " << rankNo[n] << ",  number = " << nNodes[n] << std::endl;
    }
}



#endif // LBBNDMPI_H

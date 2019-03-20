#ifndef LBBNDMPI_H
#define LBBNDMPI_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"

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
void setupBndMpi(int rank, Grid<DXQY> &grid)
{
    std::ifstream ifs_rank;
    ifs_rank.open("/home/ejette/Programs/GITHUB/badchimpp/rank.mpi");

    // READ preamble of the rank.mpi
    //  -- dimension
    std::vector<int> dim(DXQY::nD, 0);
    std::string entry_type;
    ifs_rank >> entry_type;
    for (auto &d: dim) ifs_rank >> d;
    std::cout << entry_type << " = " << dim << std::endl;

    //  -- origo
    std::vector<int> origo(DXQY::nD, 0);
    ifs_rank >> entry_type;
    for (auto &d: origo) ifs_rank >> d;
    std::cout << entry_type << " = " << origo << std::endl;

    //  -- rim
    int rim_width = 0;
    ifs_rank >> entry_type >> rim_width;
    std::cout << entry_type << " = " << rim_width << std::endl;

    std::getline(ifs_rank, entry_type);  // Read the remainer of the rim_width line
    std::getline(ifs_rank, entry_type);  // Read the <lable int> line

    std::ifstream ifs_rank_label;
    ifs_rank_label.open("/home/ejette/Programs/GITHUB/badchimpp/rank_1_labels.mpi");
    std::cout << "TEST" << std::endl;
    std::getline(ifs_rank_label, entry_type);
    std::cout << entry_type << std::endl;
    std::getline(ifs_rank_label, entry_type);
    std::cout << entry_type << std::endl;
    std::getline(ifs_rank_label, entry_type);
    std::cout << entry_type << std::endl;
    std::getline(ifs_rank_label, entry_type);
    std::cout << entry_type << std::endl;


    // Holds the node type as read from the rank.mpi file
    std::vector<int> node_type(grid.num_nodes(), 0);
    // Here we will fill the node_type with the values read from 'rank.mpi'
    //  at the node labels given by rank_#_labels.mpi


    ifs_rank.close();
    ifs_rank_label.close();

}



#endif // LBBNDMPI_H

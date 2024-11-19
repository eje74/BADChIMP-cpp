#ifndef LBMONLATMPI_H
#define LBMONLATMPI_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBboundary.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBnodes.h"
#include "../io/Input.h"
#include "../lbsolver/LButilities.h"
#include "../lbsolver/LBvtk.h"

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
 *        'Rank a node-list '[labels_ranka[r], ... ]  (Receive list)
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


class MonLatMpi
{
public:
    MonLatMpi(int neigRank, std::vector<int> &nodesToSend, std::vector<int> &nDirPerNodeToSend, std::vector<int> &dirListToSend,
              std::vector<int> &nodesReceived, std::vector<int> &nDirPerNodeReceived, std::vector<int> &dirListReceived)
        : neigRank_(neigRank), nodesToSend_(nodesToSend), nDirPerNodeToSend_(nDirPerNodeToSend), dirListToSend_(dirListToSend),
          nodesReceived_(nodesReceived), nDirPerNodeReceived_(nDirPerNodeReceived), dirListReceived_(dirListReceived) {
        // Calculate the total number of vel distirubtion to be sent to
        // the neighboring process
        std::size_t numSent = 0;
        for (auto nDir: nDirPerNodeToSend_)
            numSent += static_cast<std::size_t>(nDir);
        // if (numSent < nodesToSend_.size())
        //     numSent  = nodesToSend_.size();
        if (numSent < 3*nodesToSend_.size()) // Make sure that there is room for a 3d vector
            numSent  = 3*nodesToSend_.size();
        sendBuffer_.resize(numSent);

        std::size_t numReceived = 0;
        for (auto nDir: nDirPerNodeReceived_)
            numReceived += static_cast<std::size_t>(nDir);
        // if (numReceived < nodesReceived_.size())
        //     numReceived  = nodesReceived_.size();
        if (numReceived < 3*nodesReceived_.size()) // Make sure that there is room for a 3d vector
            numReceived  = 3*nodesReceived_.size();
        receiveBuffer_.resize(numReceived);
    }

    void inline communicateScalarField(const int &myRank, ScalarField &field, const int &fieldNo);

    template <typename DXQY>
    void inline communicateVectorField_TEST(const int &myRank, VectorField<DXQY> &field, const int &fieldNo);
    
    template <typename DXQY>
    void inline communicateLbField(const int &myRank, const Grid<DXQY> &grid, LbField<DXQY> &field, const int &fieldNo);

    void printNodesToSend() {
        std::cout << "Nodes to send to rank " << neigRank_ << ": ";
        for (auto nodeNo : nodesToSend_) {
            std::cout << " " << nodeNo;
        }
        std::cout << std::endl;
    }

    void printNodesReceived() {
        std::cout << "Nodes received from " << neigRank_ << ": ";
        for (auto nodeNo : nodesReceived_) {
            std::cout << " " << nodeNo;
        }
        std::cout << std::endl;
    }

    void printnDirToSend() {
        std::cout << "Num of dir to send to rank " << neigRank_ << " per node: ";
        for (auto nodeNo : nDirPerNodeToSend_) {
            std::cout << " " << nodeNo;
        }
        std::cout << " (" << dirListToSend_.size() << ")" << std::endl;

    }

    void printDirListToSend() {
        std::size_t cnt = 0;
        std::cout << "Dir to send to rank " << neigRank_ << std::endl;
        for (std::size_t n = 0; n < nDirPerNodeToSend_.size(); ++n) {
            for (int m = 0; m < nDirPerNodeToSend_[n]; ++m) {
                std::cout << dirListToSend_[cnt] << " ";
                cnt += 1;
            }
            std::cout << " (" << nodesToSend_[n] << ")" << std::endl;
        }

    }

    void printnDirRecived() {
        std::cout << "Num of dir recived from rank " << neigRank_ << " per node: ";
        for (auto nodeNo : nDirPerNodeReceived_) {
            std::cout << " " << nodeNo;
        }
        std::cout << " (" << dirListReceived_.size() << ")" << std::endl;
    }

    void printDirListRecived() {
        std::size_t cnt = 0;
        std::cout << "Dir recived from rank " << neigRank_ << std::endl;
        for (std::size_t n = 0; n < nDirPerNodeReceived_.size(); ++n) {
            for (int m = 0; m < nDirPerNodeReceived_[n]; ++m) {
                std::cout << dirListReceived_[cnt] << " ";
                cnt += 1;
            }
            std::cout << " (" << nodesReceived_[n] << ")" << std::endl;
        }
    }

private:
    int neigRank_; // Rank of adjacent processor
    std::vector<int> nodesToSend_;  // List of local node labels for the nodes in each mpi boundary
    std::vector<int> nDirPerNodeToSend_;  // Number of directions per node to send
    std::vector<int> dirListToSend_; // List of directions to send

    std::vector<int> nodesReceived_;  // List of which nodes the mpi boundary represents in the adjactent process
    std::vector<int> nDirPerNodeReceived_; // Number of directions per node received
    std::vector<int> dirListReceived_;  // List of direcation received

    std::vector<lbBase_t> sendBuffer_; // Buffer for sending values. sendBuffer_[i] = value(nodesToSend_[i])
    std::vector<lbBase_t> receiveBuffer_;  // Buffer for receiving values. value(nodesReceived_[i]) = receiveBuffer[i]
};


void inline MonLatMpi::communicateScalarField(const int &myRank, ScalarField &field, const int &fieldNo)
{
    if (myRank < neigRank_) {
        // SEND first
        // -- make send buffer
        for (std::size_t n=0; n < nodesToSend_.size(); ++n)
            sendBuffer_[n] = field(fieldNo, nodesToSend_[n]);
        MPI_Send(sendBuffer_.data(), static_cast<int>(nodesToSend_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD);
        // RECEIVE
        MPI_Recv(receiveBuffer_.data(), static_cast<int>(nodesReceived_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (std::size_t n=0; n < nodesReceived_.size(); ++n)
            field(fieldNo, nodesReceived_[n]) = receiveBuffer_[n];

    }  else {
        // RECEIVE first
        MPI_Recv(receiveBuffer_.data(), static_cast<int>(nodesReceived_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (std::size_t n=0; n < nodesReceived_.size(); ++n)
            field(fieldNo, nodesReceived_[n]) = receiveBuffer_[n];
        // SEND
        // -- make send buffer
        for (std::size_t n=0; n < nodesToSend_.size(); ++n)
            sendBuffer_[n] = field(fieldNo, nodesToSend_[n]);
        MPI_Send(sendBuffer_.data(), static_cast<int>(nodesToSend_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD);
    }
}

template <typename DXQY>
void inline MonLatMpi::communicateVectorField_TEST(const int &myRank, VectorField<DXQY> &vfield, const int &fieldNo)
{
  for (int d=0; d<DXQY::nD; ++d){
    if (myRank < neigRank_) {
        // SEND first
        // -- make send buffer
        for (std::size_t n=0; n < nodesToSend_.size(); ++n)
	  sendBuffer_[n] = vfield(fieldNo, d, nodesToSend_[n]);
        MPI_Send(sendBuffer_.data(), static_cast<int>(nodesToSend_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD);
        // RECEIVE
        MPI_Recv(receiveBuffer_.data(), static_cast<int>(nodesReceived_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (std::size_t n=0; n < nodesReceived_.size(); ++n)
	  vfield(fieldNo, d, nodesReceived_[n]) = receiveBuffer_[n];

    }  else {
        // RECEIVE first
        MPI_Recv(receiveBuffer_.data(), static_cast<int>(nodesReceived_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (std::size_t n=0; n < nodesReceived_.size(); ++n)
	  vfield(fieldNo, d, nodesReceived_[n]) = receiveBuffer_[n];
        // SEND
        // -- make send buffer
        for (std::size_t n=0; n < nodesToSend_.size(); ++n)
	  sendBuffer_[n] = vfield(fieldNo, d,nodesToSend_[n]);
        MPI_Send(sendBuffer_.data(), static_cast<int>(nodesToSend_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD);
    }
  }
}

template <typename DXQY>
void inline MonLatMpi::communicateLbField(const int &myRank, const Grid<DXQY> &grid, LbField<DXQY> &field, const int &fieldNo)
{

    if (myRank < neigRank_) {
        // SEND first
        // -- make send buffer
        std::size_t cnt = 0;
        for (std::size_t n=0; n < nodesToSend_.size(); ++n) {
            for (int q = 0; q < nDirPerNodeToSend_[n]; ++q) {
                int qDir = dirListToSend_[cnt];
                int ghostNode = grid.neighbor(qDir, nodesToSend_[n]);

                sendBuffer_[cnt] = field(fieldNo, qDir, ghostNode);
                cnt += 1;
            }
        }
        MPI_Send(sendBuffer_.data(), static_cast<int>(dirListToSend_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD);

        // RECEIVE

        MPI_Recv(receiveBuffer_.data(), static_cast<int>(dirListReceived_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cnt = 0;
        for (std::size_t n=0; n < nodesReceived_.size(); ++n) {
            for (int q = 0; q < nDirPerNodeReceived_[n]; ++q) {
                int qDir = dirListReceived_[cnt];
                int realNode = grid.neighbor(qDir, nodesReceived_[n]);

                field(fieldNo, qDir, realNode) = receiveBuffer_[cnt] ;
                cnt += 1;
            }
        }
    } else {
        // RECEIVE first
        MPI_Recv(receiveBuffer_.data(), static_cast<int>(dirListReceived_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::size_t cnt = 0;
        for (std::size_t n=0; n < nodesReceived_.size(); ++n) {
            for (int q = 0; q < nDirPerNodeReceived_[n]; ++q) {
                int qDir = dirListReceived_[cnt];
                int realNode = grid.neighbor(qDir, nodesReceived_[n]);

                field(fieldNo, qDir, realNode) = receiveBuffer_[cnt] ;
                cnt += 1;
            }
        }

        // SEND
        // -- Fill send buffer
        cnt = 0;
        for (std::size_t n=0; n < nodesToSend_.size(); ++n) {
            for (int q = 0; q < nDirPerNodeToSend_[n]; ++q) {
                int qDir = dirListToSend_[cnt];
                int ghostNode = grid.neighbor(qDir, nodesToSend_[n]);

                sendBuffer_[cnt] = field(fieldNo, qDir, ghostNode);
                cnt += 1;
            }
        }

        MPI_Send(sendBuffer_.data(), static_cast<int>(dirListToSend_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD);
    }
}



#endif // LBMONLATMPI_H

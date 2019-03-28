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



class MonLatMpiScalar
{
public:
    MonLatMpiScalar(int neigRank, std::vector<int> &nodesToSend, std::vector<int> &nodesReceived)
        : neigRank_(neigRank), nodesToSend_(nodesToSend), nodesReceived_(nodesReceived)
    {
        sendBuffer_.resize(nodesToSend_.size());
        receiveBuffer_.resize(nodesReceived_.size());
    }
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

private:
    int neigRank_; // List of rank of adjacent processors
    std::vector<int> nodesToSend_;  // List of local node labels for the nodes in each mpi boundary
    std::vector<int> nodesReceived_;  // List of which nodes the mpi boundary represents in the adjactent process

    std::vector<lbBase_t> sendBuffer_; // Buffer for sending values. sendBuffer_[i] = value(nodesToSend_[i])
    std::vector<lbBase_t> receiveBuffer_;  // Buffer for receiving values. value(nodesReceived_[i]) = receiveBuffer[i]
};

template <typename DXQY>
class BndMpi
{
public:
    BndMpi(int myRank): myRank_(myRank){}
    inline void addScalar(int neigRank, std::vector<int> &nodesToSend, std::vector<int> &nodesReceived)
    {
        scalarList_.emplace_back(neigRank, nodesToSend, nodesReceived);
    }
    void setupBndMpi(MpiFile<DXQY> &localFile, MpiFile<DXQY> &globalFile, MpiFile<DXQY> &rankFile, Grid<DXQY> &grid);
    void printNodesToSend();
    void printNodesRecived();
private:
    int myRank_;
    std::vector<MonLatMpiScalar> scalarList_;
};

template <typename DXQY>
void BndMpi<DXQY>::printNodesToSend()
{

    for (auto bnd : scalarList_){
        std::cout << "I'M RANK " << myRank_ <<": ";
        bnd.printNodesToSend();
        std::cout << "I'M RANK " << myRank_ <<": ";
        bnd.printNodesReceived();
    }   

}


// READ FILES
template <typename DXQY>
void BndMpi<DXQY>::setupBndMpi(MpiFile<DXQY> &localFile, MpiFile<DXQY> &globalFile, MpiFile<DXQY> &rankFile, Grid<DXQY> &grid)
/*
 * localFile  : Mpi file object that gives local labeling for the current rank
 * globalFile : Mpi file object that gives the global labeling
 * rankFile   : Mpi file object giving the ranks of the nodes. Note that a processor with
 *              rank 0 is marked with 1, since 0 is reserved as a solid marker
 * grid       : grid object for the current processor
*/
{
    std::vector<int> adjRankList; // List of rank of adjacent processors
    // std::vector<int> adjNumNodesList; // List of number of nodes in the boundary to each adjacent processor.
    std::vector<bool> nodeAdded(grid.num_nodes(), false);  // Used to register if an node has been counted

    std::vector<std::vector<int>> curProcNodeNo;  // List of local node labels where the received data, from  each mpi boundary, should be put
    std::vector<std::vector<int>> adjProcNodeNo;  // List of which nodes the mpi boundary represents in the adjactent process.


    // Find which processes that are adjacant to the current process given by 'rank'
    // And counts the number of boundary nodes in each process boundary.
    for (int pos=0; pos < static_cast<int>(localFile.size()); ++pos) {
        int nodeRank;  // Rank given in the rank.mpi file (Processor rank = nodeRank - 1)
        rankFile.getVal(nodeRank);   // Holds the rank of the node

        int localLabel, globalLabel;  // Local and global node label
        localFile.getVal(localLabel);  // Holds the local node label of the current processor
        globalFile.getVal(globalLabel);  // Holds the global labels

        if ( (nodeRank != 0) && (localLabel != 0) ) { // Not a solid and not a local default node
            nodeRank -= 1; // Get the actual rank number (A processor with rank 0 is marked with 1, since 0 is reserved as a solid marker)
            if (myRank_ != nodeRank) {  // If the rank is different from the current rank (given as the function argument 'rank')
                                     // Then it is assumend that the current node is a mpi-boundary node, and that
                                     // 'nodeRank' = rank of the adjactant processor

                // Iter : Finds the positon of the adjactant processor rank in the 'adjRankList'.
                //        If it doesn't exist iter = adjRankList.end().
                std::vector<int>::iterator iter = std::find(adjRankList.begin(), adjRankList.end(), nodeRank);

                // '..._back' means 'as the last element'
                if ( iter == adjRankList.end() ) { // If the rank is not already registered: add new rankNo and nNodes element
                   adjRankList.push_back(nodeRank);  // Adds nodeRank as the last element in adjRankList
                   // Create the procNode-lists
                   curProcNodeNo.emplace_back(1, localLabel); // Adds a std:vector<int> object with size 1 and value 'localLabel' as the last element in 'curProcNodeNo'
                   adjProcNodeNo.emplace_back(1, globalLabel);  // Adds a std:vector<int> object with size 1 and value 'globalLabel' as the last element in 'adjProcNodeNo'

                   nodeAdded[static_cast<std::size_t>(localLabel)] = true;  // Mark the local node as added
                }
                else if (!nodeAdded[static_cast<std::size_t>(localLabel)]) { // If it is not registered: increase nNodes
                    std::size_t index = static_cast<std::size_t>(std::distance(adjRankList.begin(), iter));  // Find the index of the adjacent rank
                    // Add node labels to procNode-lists
                    curProcNodeNo[index].push_back(localLabel);
                    adjProcNodeNo[index].push_back(globalLabel);

                    nodeAdded[static_cast<std::size_t>(localLabel)] = true;  // Mark the local node as added
                }
            } // End if it is an mpi boundary site
        } // End if site is Not a solid and not a local default node
    } // End for-loop reading map values


    // Add mpi boundary objects
    for (std::size_t n = 0; n < adjRankList.size(); ++n) {
        // Mission: save  'curProcNodeNo[n]' and 'adjProcNodeNo[n]' in the correct boundary objects
        // Need to send adjNumNodesList[n] to adjRankList[n]
        // need to receive adjRankList[n]'s adjNumNodesList


        // myRank_ sends adjProcNodeNo[n], that contains list of global node numbers in the adjactent process with rank adjRankList[n],
        // to rank=adjRankList[n]. These are the nodes myRank_ expect to recive data from in the future, and will be placed in myRank's
        // ghost nodes given by curProcNodeNo[n].

        if (myRank_ < adjRankList[n]) { // : send first
            // SEND
            // -- buffer size. Use tag = 0
            int bufferSizeTmp = static_cast<int>(adjProcNodeNo[n].size());
            MPI_Send( &bufferSizeTmp,  1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD);
            // -- send buffer . Use tag = 1
            MPI_Send(adjProcNodeNo[n].data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD);

             // RECEIVE
            // -- buffer size. Use tag = 0
            MPI_Recv(&bufferSizeTmp, 1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- Recive buffer. Use tag = 1
            std::vector<int> nodesToSendBuffer(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(nodesToSendBuffer.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Add the boundary to the list of mpi scalar boundaries
            addScalar(adjRankList[n], nodesToSendBuffer, curProcNodeNo[n]);


        } else if (myRank_ > adjRankList[n] ) { // : receive first
            // RECEIVE
            // -- buffer size Use tag 0
            int bufferSizeTmp;
            MPI_Recv(&bufferSizeTmp, 1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            // -- Recive buffer. Use tag = 1
            std::vector<int> nodesToSendBuffer(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(nodesToSendBuffer.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Add the boundary to the list of mpi scalar boundaries
            addScalar(adjRankList[n], nodesToSendBuffer, curProcNodeNo[n]);

            // Send buffer size
            bufferSizeTmp = static_cast<int>(adjProcNodeNo[n].size());
            MPI_Send( &bufferSizeTmp,  1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD  );
            // -- send buffer . Use tag = 1
            MPI_Send(adjProcNodeNo[n].data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD);


        } else {
            std::cout << "Error in setting up mpi boundary myRank == to adjacent rank for rank = " << myRank_ << std::endl;
            exit(EXIT_FAILURE);
        }

    } // End for all adjacent processors


    // For each Adjactant processor with rank 'adjRankList[n]'
    // Send buffer :  to the adjactant processor 'adjRankList[n]'
    //
    // INIT MPI_BOUNDARY OBJECTS
    //
    // SEND BUFFER SIZE
    // for (int n = 0; n < adjRankList.size(); ++n) {
    //   MPI_Send(
    //    &(adjProcNodeNo[n].size()),  //void* data,
    //    1, // int count,
    //    MPI_INT, // MPI_Datatype datatype,
    //    adjRankList[n], // int destination,
    //    0,  // (??)  int tag,
    //    MPI_COMM_WORLD  //(??) MPI_Comm communicator
    //   )
    // }
    //
    // RECEIVE BUFFER SIZE
    //
    // for (int n = 0; n < adjRankList.size(); ++n) {
    //   int bufferSize;
    //   MPI_Recv(
    //    &bufferSize //void* data,
    //    1, //int count,
    //    MPI_INT, // MPI_Datatype datatype,
    //    adjRankList[n], // int source,
    //    0, // (??) int tag,
    //    MPI_COMM_WORLD,  // MPI_Comm communicator,
    //    MPI_STATUS_IGNORE  //  MPI_Status* status
    //   )
    //   -- setup recive boundary size and buffer
    // }

    //
    // for (int n = 0; n < adjRankList.size(); ++n) {
    //   MPI_Send(
    //    adjProcNodeNo[n],  //void* data,
    //    adjProcNodeNo[n].size(), // int count,
    //    MPI_INT, // MPI_Datatype datatype,
    //    adjRankList[n], // int destination,
    //    1,  //  int tag,
    //    MPI_COMM_WORLD (??) //MPI_Comm communicator
    //   )
    // }

    //
    //
    // for (int n = 0; n < adjRankList.size(); ++n) {
    //   MPI_Recv(
    //    &buffer, //void* data,
    //    bufferSize[n], //int count,
    //    MPI_INT, // MPI_Datatype datatype,
    //    adjRankList[n],
    //    1, // int tag,
    //    MPI_COMM_WORLD,  // MPI_Comm communicator,
    //    MPI_STATUS_IGNORE  //  MPI_Status* status
    //   )
    // }
    //
    // How to fill a scalar mpi send buffer?
    // sendRhoBuffer[m] = rho(0, buffer[m])
    //
    // send to adj proc
    // rho(0, curProcNodeNo[m]) = from adj's sendRhoBuffer[m]
}



#endif // LBBNDMPI_H

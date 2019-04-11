#ifndef LBBNDMPI_H
#define LBBNDMPI_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
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
          nodesReceived_(nodesReceived), nDirPerNodeReceived_(nDirPerNodeReceived), dirListReceived_(dirListReceived)
    {
        // Calculate the total number of vel distirubtion to be sent to
        // the neighboring process
        std::size_t numSent = 0;
        for (auto nDir: nDirPerNodeToSend)
            numSent += static_cast<std::size_t>(nDir);
        if (numSent < nodesToSend_.size())
            numSent  = nodesToSend_.size();
        sendBuffer_.resize(numSent);

        std::size_t numReceived = 0;
        for (auto nDir: nDirPerNodeToSend)
            numReceived += static_cast<std::size_t>(nDir);
        if (numReceived < nodesToSend_.size())
            numReceived  = nodesToSend_.size();
        receiveBuffer_.resize(numReceived);
    }

    void inline communicateScalarField(const int &myRank, ScalarField &field);

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


void inline MonLatMpi::communicateScalarField(const int &myRank, ScalarField &field)
{
    if (myRank < neigRank_) {
        // SEND first
        // -- make send buffer
        for (std::size_t n=0; n < nodesToSend_.size(); ++n)
            sendBuffer_[n] = field(0, nodesToSend_[n]);
        MPI_Send(sendBuffer_.data(), static_cast<int>(nodesToSend_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD);

        MPI_Recv(receiveBuffer_.data(), static_cast<int>(nodesReceived_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (std::size_t n=0; n < nodesReceived_.size(); ++n)
            field(0, nodesReceived_[n]) = receiveBuffer_[n];

    }  else {
        // RECEIVE first
        MPI_Recv(receiveBuffer_.data(), static_cast<int>(nodesReceived_.size()), MPI_DOUBLE, neigRank_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (std::size_t n=0; n < nodesReceived_.size(); ++n)
            field(0, nodesReceived_[n]) = receiveBuffer_[n];

        // -- make send buffer
        for (std::size_t n=0; n < nodesToSend_.size(); ++n)
            sendBuffer_[n] = field(0, nodesToSend_[n]);
        MPI_Send(sendBuffer_.data(), static_cast<int>(nodesToSend_.size()), MPI_DOUBLE, neigRank_, 1, MPI_COMM_WORLD);
    }
}


template <typename DXQY>
void inline MonLatMpi::communicateLbField(const int &myRank, const Grid<DXQY> &grid, LbField<DXQY> &field, const int &fieldNo) {

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
                int qDir = dirListToSend_[cnt];
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
                int qDir = dirListToSend_[cnt];
                int realNode = grid.neighbor(qDir, nodesReceived_[n]);
                field(fieldNo, qDir, realNode) = receiveBuffer_[cnt] ;
                cnt += 1;
            }
        }
        // SEND
        // -- make send buffer
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



template <typename DXQY>
class BndMpi
{
public:
    BndMpi(int myRank): myRank_(myRank){}

    inline void addMpiBnd(int neigRank, std::vector<int> &nodesToSend, std::vector<int> &nDirPerNodeToSend, std::vector<int> &dirListToSend,
                    std::vector<int> &nodesReceived, std::vector<int> &nDirPerNodeReceived, std::vector<int> &dirListReceived)
    {
        //scalarList_.emplace_back(neigRank, nodesToSend, nodesReceived);
        mpiList_.emplace_back(neigRank, nodesToSend, nDirPerNodeToSend, dirListToSend, nodesReceived, nDirPerNodeReceived, dirListReceived);
    }
    void inline communciateScalarField(ScalarField &field);
    void inline communicateLbField(Grid<DXQY> &grid, LbField<DXQY> &field, const int fieldNo);
    void setupBndMpi(MpiFile<DXQY> &localFile, MpiFile<DXQY> &globalFile, MpiFile<DXQY> &rankFile, Grid<DXQY> &grid);
    void printNodesToSend();
    void printNodesRecived();
private:
    int myRank_;
    std::vector<MonLatMpi> mpiList_;
};

template <typename DXQY>
void inline BndMpi<DXQY>::communciateScalarField(ScalarField &field)
{
    for (auto& mpibnd: mpiList_)
        mpibnd.communicateScalarField(myRank_, field);
}

template <typename DXQY>
void inline BndMpi<DXQY>::communicateLbField(Grid<DXQY> &grid, LbField<DXQY> &field, const int fieldNo)
{
    for (auto& mpibnd: mpiList_)
        mpibnd.communicateLbField(myRank_, grid, field, fieldNo);
}


template <typename DXQY>
void BndMpi<DXQY>::printNodesToSend()
{

    for (auto bnd : mpiList_){
    //for (auto bnd : scalarList_){
        std::cout << "I'M RANK " << myRank_ <<": ";
        bnd.printNodesToSend();
        std::cout << "I'M RANK " << myRank_ <<": ";
        bnd.printNodesReceived();
        std::cout << "I'M RANK " << myRank_ <<": ";
        bnd.printnDirToSend();
        std::cout << "I'M RANK " << myRank_ <<": ";
        bnd.printDirListToSend();
    }

}

template <typename DXQY>
void makeDirList(int myRank, Grid<DXQY> &grid,  std::vector<int> &curProcNodeNo,  std::vector<int> &nDirPerNode, std::vector<int> &dirList)
{
    // set size of the number of directino per node vector. Same as the number of ghost nodes
    nDirPerNode.resize(curProcNodeNo.size());
    for (std::size_t n = 0; n < curProcNodeNo.size(); ++n) { // List of ghost nodes between two adjaction processors
        int currentGhostNode = curProcNodeNo[n];
        int nDir = 0;  // Counter for number of directions per ghost node


        for (int q = 0; q < DXQY::nQNonZero_; ++q) { // Loop over all directions exept the rest direction
            int ghostNodeNeig = grid.neighbor(q, currentGhostNode);
            if ( grid.getRank(ghostNodeNeig) == myRank) {  // Direction from ghost node to bulk node of the current rank
                // This value is part of the mpi boundary
                dirList.push_back(q);
                nDir += 1;
            }
        }  // END Loop through all directions
        nDirPerNode[n] = nDir;
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
    std::vector<bool> nodeAdded(grid.size(), false);  // Used to register if an node has been counted

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

            // SEND microscopic fileds mpi info
            std::vector<int> nDirPerNode;  // List of directions we want to recive from the adjactent node
            std::vector<int> dirList;  // List of directions we want do recive
            makeDirList(myRank_, grid,  curProcNodeNo[n], nDirPerNode, dirList);

            // -- nDirPerNode. use tag = 2
            MPI_Send(nDirPerNode.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 2, MPI_COMM_WORLD);
            // -- dirList.size use tag = 3
            bufferSizeTmp = static_cast<int>(dirList.size());
            MPI_Send( &bufferSizeTmp,  1, MPI_INT, adjRankList[n], 3, MPI_COMM_WORLD);
            // -- dirList use tag = 4
            MPI_Send(dirList.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 4, MPI_COMM_WORLD);

             // RECEIVE
            // -- buffer size. Use tag = 0
            MPI_Recv(&bufferSizeTmp, 1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- Recive buffer. Use tag = 1
            std::vector<int> nodesToSendBuffer(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(nodesToSendBuffer.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- nDirPerNode. use tag = 2
            std::vector<int> nDirPerNodeToSend(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(nDirPerNodeToSend.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- dirList.size use tag = 3
            MPI_Recv(&bufferSizeTmp, 1, MPI_INT, adjRankList[n], 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- dirList use tag = 4
            std::vector<int> dirListToSend(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(dirListToSend.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Add the boundary to the list of mpi scalar boundaries
            addMpiBnd(adjRankList[n], nodesToSendBuffer, nDirPerNodeToSend, dirListToSend, curProcNodeNo[n], nDirPerNode, dirList);


        } else if (myRank_ > adjRankList[n] ) { // : receive first
            // RECEIVE
            // -- buffer size Use tag 0
            int bufferSizeTmp;
            MPI_Recv(&bufferSizeTmp, 1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            // -- Recive buffer. Use tag = 1
            std::vector<int> nodesToSendBuffer(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(nodesToSendBuffer.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- nDirPerNode. use tag = 2
            std::vector<int> nDirPerNodeToSend(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(nDirPerNodeToSend.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- dirList.size use tag = 3
            MPI_Recv(&bufferSizeTmp, 1, MPI_INT, adjRankList[n], 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // -- dirList use tag = 4
            std::vector<int> dirListToSend(static_cast<std::size_t>(bufferSizeTmp));
            MPI_Recv(dirListToSend.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Make the nDirPerNode and dirList
            std::vector<int> nDirPerNode;  // List of directions we want to recive from the adjactent node
            std::vector<int> dirList;  // List of directions we want do receive
            makeDirList(myRank_, grid, curProcNodeNo[n], nDirPerNode, dirList);

            // Add the boundary to the list of mpi scalar boundaries
            addMpiBnd(adjRankList[n], nodesToSendBuffer, nDirPerNodeToSend, dirListToSend, curProcNodeNo[n], nDirPerNode, dirList);
            // SEND
            // --  Send buffer size. Use tag = 0
            bufferSizeTmp = static_cast<int>(adjProcNodeNo[n].size());
            MPI_Send( &bufferSizeTmp,  1, MPI_INT, adjRankList[n], 0, MPI_COMM_WORLD  );
            // -- send buffer . Use tag = 1
            MPI_Send(adjProcNodeNo[n].data(), bufferSizeTmp, MPI_INT, adjRankList[n], 1, MPI_COMM_WORLD);
            // SEND microscopic fields mpi info
            // -- nDirPerNode. use tag = 2
            MPI_Send(nDirPerNode.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 2, MPI_COMM_WORLD);
            // -- dirList.size use tag = 3
            bufferSizeTmp = static_cast<int>(dirList.size());
            MPI_Send( &bufferSizeTmp,  1, MPI_INT, adjRankList[n], 3, MPI_COMM_WORLD);
            // -- dirList use tag = 4
            MPI_Send(dirList.data(), bufferSizeTmp, MPI_INT, adjRankList[n], 4, MPI_COMM_WORLD);


        } else {
            std::cout << "Error in setting up mpi boundary myRank == to adjacent rank for rank = " << myRank_ << std::endl;
            exit(EXIT_FAILURE);
        }

    } // End for all adjacent processors

    // Reset file for reading
    localFile.reset();
    globalFile.reset();
    rankFile.reset();
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

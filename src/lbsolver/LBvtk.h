/*
 * LVvtk.h
 *
 *  Created on 14. September
 *      Author: esje
 */

#ifndef SRC_LBVTK_H_
#define SRC_LBVTK_H_

#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <typeinfo>
#include <numeric>
#include <algorithm>
#include <set>
#include <map>
#include "LBlatticetypes.h"
//#include "LButilities.h"


/* *********************** TO DO ***********************
 * - add number of components to the data set attributes
 * - add point data sub set key word POINT_DATA_SUBSET
 *
 */


/* A MpiFile object is used to read the files generate by the
   mpiGrid.py script.
   it will hold the premable information
   dim_      : dimensions including the rim
   rimWidth_ : with of the rim (in number of nodes)
   origo_    : position of the local origo (so that the global pos is 'local pos' + 'origo_'

   It also contains the functions
   getVal(var)  : reads one label/rank value from the map
   size()       : size of the label/rank map
   reset()      : resets the file head so that it will read from
                  the beginning of the label/rank map
*/
template<typename DXQY>
class LBvtk
{
public:
    LBvtk(std::string filename);

    ~LBvtk() {
        ifs_.close();
    }

    void readPreamble();
    std::string readDataSetType();
    void readUNSTRUCTURED_LB_GRID();
    void readPoints();
    void toPos() {
        ifs_.seekg(beginPosBlock_, std::ios_base::beg);
    }
    template<typename T>
    std::vector<T> getPos();

    void readLattice();
    void readNeighbors();
    void toNeighbors() {
        ifs_.seekg(beginNeigBlock_, std::ios_base::beg);
    }
    template<typename T>
    std::vector<T> getNeighbors();

    void readMpi();
    void readPointData();
    void toAttribute(const std::string &dataName);
    template<typename T>
    T getScalarAttribute() {
        T ret;
        ifs_ >> ret;
        return ret;
    }
    template<typename T>
    T getScalar() {
        T ret;
        ifs_ >> ret;
        return ret;
    }

/*    template<typename T>
    struct dataSubsetElement
    {
        T val;
        int nodeNo;        
    };
    

    template<typename T>
    struct dataSubsetElement<T> getSubsetAttribure() {
        struct dataSubsetElement<T> ret;
        ifs_ >> ret.nodeNo >> ret.val;
        return ret;
    } */
    

    template<typename T>
    auto  getSubsetAttribure() {
        struct dataSubsetElement
        {
            T val;
            int nodeNo;        
        } ret;
        ifs_ >> ret.nodeNo >> ret.val;
        return ret;
    }


    void readPointDataSubset();

    inline int numSubsetEntries() const {
        return numSubsetEntries_;
    }

    inline int getNumPoints() const {
        return nPoints_;
    }


    inline int getGlobaDimensions(const int d) const {
        return globalDimensions_[d];
    }
    
    inline std::vector<int> getGlobaDimensions() const {
        std::vector<int> dim(3,1);
        // std::vector<int> dim(DXQY::nD,1);
        for (int i=0; i < DXQY::nD; ++i)
            dim[i] = globalDimensions_[i];
        return dim;
    }

    inline int beginNodeNo() const {
        return zero_ghost_node_ ? 1 : 0;
    }

    inline int endNodeNo() const {
        return zero_ghost_node_ ? (nPoints_ + 1) : nPoints_;
    }

    // Number of neighboring processes
    int getNumNeigProc() const {
        return adjRankList_.size();
    }

    inline int getRank() const {return rank_;}

    // Get the rank of neighboring processor n
    int getNeigRank(int const n) const {
        return adjRankList_[n];
    }

    // Get the node numbers on the current processor of
    // the nodes that overlapp with neighbor processor n
    std::vector<int> getNeigNodesNo(const int n) const {
        return curProcNodeNo_[n];
    }

    // Get the node numbers on neighbor processor n of
    // the nodes that overlapp with this one
    std::vector<int> getNeigNodesNeigNo(const int n) const {
        return adjProcNodeNo_[n];
    }

    inline int dirVtkToLB(const int q) {return f2p_[q];}

private:
    std::string filename_;  // The file name
    std::ifstream ifs_;   // File stream object

    // PREAMBLE
    std::string title_; // File title
    bool ascii_; // Data type for file (binary or ascii)

    // DATASET
    std::set<std::string> dataTypeList_ {"UNSTRUCTURED_LB_GRID"};
    //std::set<std::string> dataSetAttribute_ {"SCALARS"}
    bool zero_ghost_node_;
    std::vector<int> globalDimensions_;
    int nD_; // Number of spatial dimensions

    // POINTS
    int nPoints_; // Number of points
    int beginPosBlock_;
    std::string dataType_; // Current datatype

    // LATTICE
    int nQ_; // Number of basis vector
    std::vector<int> f2p_; // File to program: maps the basis
    // used in the input file to basis used in the program
    int beginNeigBlock_;

    // USED for initiating mpi transfers
    int rank_; // Rank of the process
    std::vector<std::vector<int>> curProcNodeNo_;  // List of local node labels where the received data, from  each mpi boundary, should be put
    std::vector<std::vector<int>> adjProcNodeNo_;  // List of which nodes the mpi boundary represents in the adjactent process.
    std::vector<int> adjRankList_;  // List of ranks of the neighboring processes

    // DATA SET ATTRIBUTES
    std::map<std::string, int> dataAttributes_;

    // DATA SUBSET ATTRUBUTES
    std::map<std::string, std::vector<long int>> dataSubsetAttributes_; // block name, [size, file position]

    int numSubsetEntries_;
    enum { UNDEFINED, POINT_DATA, POINT_DATA_SUBSET } readState_;

};


template<typename DXQY>
LBvtk<DXQY>::LBvtk(std::string filename) : filename_(filename)
{
    // Set default values
    nD_ = 3; // Default value
    zero_ghost_node_ = false; // Default value

    // Open input file
    ifs_.open(filename_);
    if (!ifs_) {
        std::cout << "ERROR reading file in LBvtk: coult not open file " << filename_ << "." << std::endl;
        exit(1);
    }

    // Preamble
    readPreamble();

    // Read data set
    std::string dataSetType = readDataSetType();
    if (dataSetType == "UNSTRUCTURED_LB_GRID") {
        readUNSTRUCTURED_LB_GRID();
        // POINTS
        readPoints();
        // LATTICE
        readLattice();
        // NEIGHBORS
        readNeighbors();
    }

    // Parallell computing
    // READ MPI
    readMpi();

    // Read data attribute
    readPointData();

    // Read data subset attribues
    readPointDataSubset();

    readState_ = UNDEFINED;
    // Since the file is read to the end we need to reopen the file.
}


template<typename DXQY>
void LBvtk<DXQY>::readPreamble()
{
    std::string str;

    // Header
    std::getline(ifs_, str);

    // Titel
    ifs_ >> title_;
    std::getline(ifs_, str);

    // ASCII | BINARY  <data type>
    ifs_ >> str;
    // -- set datatype
    if (str == "ASCII") {
        ascii_ = true;
    } else if  (str == "BINARY") {
        ascii_ = false;
    } else { // throw error.
        std::cout << "Files data type must be ASCII or BINARY" << std::endl;
        exit(1);
    }
    // Read the rest of the line
    std::getline(ifs_, str);
}


template<typename DXQY>
std::string LBvtk<DXQY>::readDataSetType()
{
    std::string str, tmp;

    // DATASET type   <geometry/topology>
    // Read keyword
    ifs_ >> str;
    if ( str != "DATASET" ) {
        std::cout << "Expected keyword: DATA_TYPE" << std::endl;
        exit(1);
    }
    // Read data set type
    ifs_ >> str;
    if (dataTypeList_.find(str) == dataTypeList_.end()) {
        std::cout << "Unknown DATASET type: " << str << std::endl;
        exit(1);
    }
    // Read rest of line
    std::getline(ifs_, tmp);

    return str;
}


template<typename DXQY>
void LBvtk<DXQY>::readUNSTRUCTURED_LB_GRID()
{
    std::string str;

    // Two optional keywords
    // NUM_DIMENSIONS nd  <nd: int number of spatial dimension> [OPTIONAL]
    int ifs_pos = ifs_.tellg();
    ifs_ >> str;
    if ( str == "NUM_DIMENSIONS") {
        ifs_ >> nD_;
        std::getline(ifs_, str);
    } else {
        // Reset file
        ifs_.seekg(ifs_pos ,std::ios_base::beg);
    }

    ifs_ >> str;
    if ( str == "GLOBAL_DIMENSIONS" ) {
        for (int d = 0; d < nD_; ++d) {
            int val;
            ifs_ >> val;
            globalDimensions_.push_back(val);
        }
        std::getline(ifs_, str);
    } else {
        std::cout << "Expected keyword: DATA_TYPE got " << str << std::endl;
        exit(1);
    }


    // USE_ZERO_GHOST_NODE  <use node 0 as the default ghost node> [OPTIONAL]
    ifs_pos = ifs_.tellg();
    ifs_ >> str;
    if ( str == "USE_ZERO_GHOST_NODE" ) {
        zero_ghost_node_ = true;
        std::getline(ifs_, str);
    } else {
        // Reset file
        ifs_.seekg(ifs_pos ,std::ios_base::beg);
    }
}


template<typename DXQY>
void LBvtk<DXQY>::readPoints()
{
    std::string str;
    // POINTS n dataType <n: number of points>
    ifs_ >> str;
    if ( str == "POINTS" ) {
        ifs_ >> nPoints_ >> dataType_;
        std::getline(ifs_, str); // Read end of line
    } else {
        std::cout << "Unexptected keyword: " << str << std::endl;
        exit(1);
    }
    beginPosBlock_ = ifs_.tellg();

    // Read through the rest of the points
    for (int n = 0; n < (nD_*nPoints_); ++n) {
        int tmp;
        ifs_ >> tmp;
    }
    std::getline(ifs_, str);
}


template<typename DXQY>
template<typename T>
std::vector<T> LBvtk<DXQY>::getPos() 
{
    std::vector<T> pos(nD_, 0);
    for (auto &x : pos) {
        ifs_ >> x;
    }
    return pos;
}


template<typename DXQY>
void LBvtk<DXQY>::readLattice()
{
    std::string str;

    ifs_ >> str;

    if ( str == "LATTICE") {
        ifs_ >> nQ_ >> dataType_;
        std::getline(ifs_, str);
    } else {
        std::cout << "Unexptected entry. Expected key word LATTICE got  " << str << std::endl;
        exit(1);
    }

    if (nQ_ != DXQY::nQ) {
        std::cout << "Error reading BADChIMP vtklb input file: " << filename_ << ", in section LATTICE.\n Wrong number of basis vectors. Expected " << DXQY::nQ << " got " << nQ_ << std::endl;
        exit(1);
    }

    f2p_.reserve(nQ_); // Map file basis to program/header basis
    std::vector<bool> basis_reg(nQ_, false); // Check that basis vector in file is one to one with basis in header

    for ( int q = 0; q < nQ_; ++q ) {
        std::vector<int>  vec(nD_); // Basis vector
        for (auto &x : vec) {
            ifs_ >> x;
        }

        // Map file basis to program basis
        f2p_[q] = DXQY::c2q(vec);
        basis_reg[f2p_[q]] = true;
    }

    // All basis_vec should be true
    for (auto x : basis_reg) {
        if (!x) {
            std::cout << "The LATTICE basis in file " << filename_ << " is not complete" << std::endl;
            exit(1);
        }
    }
    std::getline(ifs_, str);
}


template<typename DXQY>
void LBvtk<DXQY>::readNeighbors()
{
    std::string str;
    ifs_>> str;

    // NEIGHBORS dataType <list of neighbor nodes for each node>
    if ( str == "NEIGHBORS") {
        ifs_ >> dataType_;
        std::getline(ifs_, str);
    } else {
        std::cout << "Unrecognized keyword. Expected NEIGHBORS got " << str << std::endl;
        exit(1);
    }

    // Beginning of the neighbors block
    beginNeigBlock_ = ifs_.tellg();

    // Read
    for (int n = 0; n < (nQ_*nPoints_); ++n) {
        int tmp;
        ifs_ >> tmp;
    }

    std::getline(ifs_, str);
}


template<typename DXQY>
template<typename T>
std::vector<T> LBvtk<DXQY>::getNeighbors() 
{
    std::vector<T> neig(nQ_, 0);
    for (int q = 0; q < nQ_; ++q) {
        T x;
        ifs_ >> x;
        neig[f2p_[q]] = x;
    }
    return neig;
}


template<typename DXQY>
void LBvtk<DXQY>::readMpi()
{
    // PARALLEL_COMPUTING rank  <rank: rank of the current processor>
    std::string str;

    ifs_ >> str;
    if ( str == "PARALLEL_COMPUTING" ) {
        ifs_ >> rank_;
        std::getline(ifs_, str);
    } else {
        std::cout << "Unexptected entry. Expected PARALLEL_COMPUTING, but got " << str << std::endl;
        exit(1);
    }

    // CAN have borders with multiple processors
    // READ ALL PROCESSOR BLOCKS
    while (1) {
        // PROCESSOR n rank  <n: number of point, rank: neighbor processor>
        int numNeig;
        int rankNeig;
        // Peek at the key word
        // Get the current file position
        int ifs_pos = ifs_.tellg();
        ifs_ >> str;
        if ( str == "PROCESSOR" ) {
            ifs_ >> numNeig >> rankNeig;
            std::getline(ifs_, str);
        } else {
            // Reset file position
            ifs_.seekg(ifs_pos ,std::ios_base::beg);
            break; // Break out of the while--loop
        }

        // Read the neighbor pairs
        std::vector<int> curNodes(numNeig), adjNodes(numNeig);
        for ( int n = 0; n < numNeig; ++n ) {
            int nodeNo, neigNo;
            ifs_ >> nodeNo >> neigNo;
            curNodes[n] = nodeNo;
            adjNodes[n] = neigNo;
        }
        adjRankList_.push_back(rankNeig);
        curProcNodeNo_.push_back(curNodes);
        adjProcNodeNo_.push_back(adjNodes);

        // Read rest of line
        std::getline(ifs_, str);
    } // END of while loop block reading
    // Sort the neighborlists
    // Here we need to sort the lists
    auto idx = sort_indexes(adjRankList_);
    auto adjRankListCopy = adjRankList_;
    auto curProcNodeNoCopy = curProcNodeNo_;
    auto adjProcNodeNoCopy = adjProcNodeNo_;
    for (unsigned i = 0; i < idx.size(); ++i) {
        adjRankList_[i] = adjRankListCopy[idx[i]];
        curProcNodeNo_[i] = curProcNodeNoCopy[idx[i]];
        adjProcNodeNo_[i] = adjProcNodeNoCopy[idx[i]];
    }
}


template<typename DXQY>
void LBvtk<DXQY>::readPointData()
{
    std::string str;
    int getNumPoints;

    // POINT_DATA n    <dataset attributes>
    ifs_ >> str;

    if (str == "POINT_DATA") {
        ifs_ >> getNumPoints;
    } else {
        std::cout << "Unrecognized keyword. Expected POINT_DATA got " << str << std::endl;
        exit(1);
    }
    std::getline(ifs_, str); // Read resf of line

    // Read the data set attributes
    int ifs_pos = ifs_.tellg();
    while(ifs_ >> str) {         // Peak at value
        if (str == "SCALARS") {    // SCALARS dataName dataType
            std::string dataName, dataType;
            ifs_ >> dataName >> dataType;

            std::getline(ifs_, str); // Read rest of line

            // Add attribute to map
            dataAttributes_.insert(std::pair<std::string, int>(dataName, ifs_.tellg()));

            // Read the scalar entries
            for (int n=0; n < getNumPoints; ++n) getScalarAttribute<double>();
            std::getline(ifs_, str); // Read rest of line
        } else {
            // Reset file position
            ifs_.seekg(ifs_pos ,std::ios_base::beg);
            break; // Break out of the while--loop
        }
        ifs_pos = ifs_.tellg();
    }
}


template<typename DXQY>
void LBvtk<DXQY>::toAttribute(const std::string &dataName)
{
    // Check if entry exists
    std::map<std::string, int>::iterator it;
    it = dataAttributes_.find(dataName);
    if (it != dataAttributes_.end()) {
        ifs_.seekg(dataAttributes_[dataName], std::ios_base::beg);
        readState_ = POINT_DATA;
        return;
    }
    std::map<std::string, std::vector<long int>>::iterator it2;
    it2 = dataSubsetAttributes_.find(dataName);
    if (it2 != dataSubsetAttributes_.end()) {
        const auto val = dataSubsetAttributes_[dataName];
        readState_ = POINT_DATA_SUBSET;
        numSubsetEntries_ = val[0];
        const long int filePos = val[1];
        ifs_.seekg(filePos, std::ios_base::beg);
        return;
    }
    std::cout << "Data attribute with name " << dataName << " do not exist." << std::endl;
    exit(1);
    return;
}


template<typename DXQY>
void LBvtk<DXQY>::readPointDataSubset()
{
    std::string str;
    // -------------------------------------------------------------------------- ADD CODE END
    if (ifs_.eof()) { // Reset file 'end of file' so we can read from it later
        ifs_.clear();
    } 
    else {
        // POINT_DATA    <dataset attributes>
        ifs_ >> str;
        if (str == "POINT_DATA_SUBSET") {
            std::getline(ifs_, str); // Read resf of line
            // READ    
            // Read the data set attributes
            int ifs_pos = ifs_.tellg();
            while(ifs_ >> str) {         // Peak at value
                if (str == "SCALARS") {    // SCALARS dataName dataType
                    int size;
                    std::string dataName, dataType;
                    ifs_ >> size >> dataName >> dataType;
                    std::getline(ifs_, str); // Read rest of line

                    // Add attribute to map
                    dataSubsetAttributes_.insert( std::pair< std::string, std::vector<long int> >(dataName, std::vector<long int>{size, ifs_.tellg()}) );

                    // Read the scalar entries
                    for (int n=0; n < size; ++n) {
                        std::getline(ifs_, str); // Read rest of line
                    }
                } else {
                    // Reset file position
                    ifs_.seekg(ifs_pos ,std::ios_base::beg);
                    break; // Break out of the while--loop
                }
                ifs_pos = ifs_.tellg();
            }
        }

        if (ifs_.eof()) { // Reset file 'end of file' so we can read from it later
            ifs_.clear();
        }
    }
}


#endif /* SRC_LBVTK_H_ */

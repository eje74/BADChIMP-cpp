#ifndef LBBOUNDARYNEW_H
#define LBBOUNDARYNEW_H

#include <cassert>
#include <vector>
#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBnodes.h"
#include "LBgrid.h"



//=====================================================================================
//                            B O U N D A R Y   N O D E
//=====================================================================================
//  - helper class for the class Boundary 
//-------------------------------------------------------------------------------------
template<typename DXQY>
class BoundaryNode
{
public:
    BoundaryNode(){}
    template<typename T>
    BoundaryNode(const int nodeNo, const T& betaList, const T& gammaList, const T& deltaList)
    :nodeNo_(nodeNo), nBeta_(betaList.size()), nGamma_(gammaList.size()), nDelta_(deltaList.size()), 
     nUnknown_(nBeta_+2*nDelta_), deltaBegin_(0), betaBegin_(2*nDelta_), 
     gammaBegin_(betaBegin_ + 2*nBeta_)
    {
        int cnt = 0;

        auto append_to_list = [&cnt, this](const T& input_list) {
            for (const auto &dir: input_list) 
                this->nodeList_[cnt++] = dir;
            for (const auto &dir: input_list) 
                this->nodeList_[cnt++] = DXQY::reverseDirection(dir);
        };

        append_to_list(deltaList);
        append_to_list(betaList);
        append_to_list(gammaList);
    }

    inline int nodeNo() const {return nodeNo_;}
    inline int nBeta() const {return nBeta_;}
    inline int nGamma() const {return nGamma_;}
    inline int nDelta() const {return nDelta_;}
    inline std::vector<int> beta() const {return std::vector(nodeList_ + betaBegin_, nodeList_ + betaBegin_ + nBeta_);}
    inline std::vector<int> betaRev() const {return std::vector(nodeList_ + betaBegin_ + nBeta_, nodeList_ + gammaBegin_);}
    inline std::vector<int> gamma() const {return std::vector(nodeList_ + gammaBegin_, nodeList_ + gammaBegin_ + nGamma_);}
    inline std::vector<int> gammaRev() const {return std::vector(nodeList_ + gammaBegin_ + nGamma_, nodeList_ + DXQY::nQ-1);}
    inline std::vector<int> delta() const {return std::vector(nodeList_, nodeList_ + nDelta_);}        
    inline std::vector<int> deltaRev() const {return std::vector(nodeList_ + nDelta_, nodeList_ + betaBegin_);}    
    inline std::vector<int> known() const {return std::vector(nodeList_ + nUnknown_, nodeList_ + DXQY::nQ-1);}    
    inline std::vector<int> unknown() const {return std::vector(nodeList_, nodeList_ + nUnknown_);}    

private:
    const int nodeNo_;
    const int nBeta_;
    const int nGamma_;
    const int nDelta_;
    const int nUnknown_;
    const int deltaBegin_;
    const int betaBegin_;
    const int gammaBegin_;
    int nodeList_[DXQY::nQ-1];
};

//=====================================================================================
//
//                                 B O U N D A R Y
//
//=====================================================================================
template <typename DXQY>
class Boundary
{
public:
    Boundary() {}
    Boundary(const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);

    inline std::vector<int> beta(const int bndNo) const {return boundaryNodes_[bndNo].beta();}
    inline std::vector<int> gamma(const int bndNo) const {return boundaryNodes_[bndNo].gamma();}
    inline std::vector<int> delta(const int bndNo) const {return boundaryNodes_[bndNo].delta();}
    inline std::vector<int> betaRev(const int bndNo) const {return boundaryNodes_[bndNo].betaRev();}
    inline std::vector<int> gammaRev(const int bndNo) const {return boundaryNodes_[bndNo].gammaRev();}
    inline std::vector<int> deltaRev(const int bndNo) const {return boundaryNodes_[bndNo].deltaRev();}
    inline std::vector<int> unknown(const int bndNo) const {return boundaryNodes_[bndNo].unknown();}
    inline std::vector<int> known(const int bndNo) const {return boundaryNodes_[bndNo].known();}
    inline int dirRev(const int dir) const;
    inline int nodeNo(const int bndNo) const {return boundaryNodes_[bndNo].nodeNo();}
    inline int size() const {return boundaryNodes_.size();}
    inline int nBeta(const int bndNo) const {return boundaryNodes_[bndNo].nBeta();}
    inline int nGamma(const int bndNo) const {return boundaryNodes_[bndNo].nGamma();}
    inline int nDelta(const int bndNo) const {return boundaryNodes_[bndNo].nDelta();}
    
    const std::vector<BoundaryNode<DXQY>> & operator () () const {return boundaryNodes_;}
    const BoundaryNode<DXQY> & operator () (int bndNo) const {return boundaryNodes_[bndNo];}

protected:
    std::vector<BoundaryNode<DXQY>> boundaryNodes_;

    auto getLatticePairs(std::vector<bool> const & isSolid) const;

};


//                                  Boundary
//------------------------------------------------------------------------------------- getLatticePairs
template<typename DXQY>
auto Boundary<DXQY>::getLatticePairs(std::vector<bool> const & isSolid) const
{
    struct {
        std::vector<int> delta;
        std::vector<int> beta;
        std::vector<int> gamma;
    } ret;

    std::vector<bool> added(DXQY::nQ-1, false);
    for (int q = 0; q < DXQY::nQ-1; ++q) {
        if (!added[q]) {
            const int qRev = DXQY::reverseDirection(q);
            added[q] = true;
            added[qRev] = true;
            int numSolids = 0;
            numSolids += isSolid[q] ? 1 : 0; 
            numSolids += isSolid[qRev] ? 1 : 0; 

            if (numSolids == 0) {
                ret.gamma.push_back(q);
            } 
            else if (numSolids == 1) {
                const int qAdd = isSolid[qRev] ? q : qRev;
                ret.beta.push_back(qAdd);
            }
            else if (numSolids == 2) {
                ret.delta.push_back(q);
            } 
            else {
                std::cout << "Error in Boundary: number of solids in a pair is either 0, 1, 2 but is now " << numSolids << '\n';
                exit(1);
            }
        }
    }
    return ret;
}


//                                  Boundary
//------------------------------------------------------------------------------------- Boundary
template <typename DXQY>
Boundary<DXQY>::Boundary(const std::vector<int> &bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
{
    for (const auto & nodeNo: bndNodes) {
        std::vector<bool> isSolid(DXQY::nQ-1);
        for (int q = 0; q < DXQY::nQ-1; ++q) {
            const int neigNo = grid.neighbor(q, nodeNo);
            isSolid[q] = nodes.isSolid(neigNo);
        }
        auto latPrs = getLatticePairs(isSolid);
        boundaryNodes_.emplace_back(nodeNo, latPrs.beta, latPrs.gamma, latPrs.delta);
    }
}

//                                  Boundary
//------------------------------------------------------------------------------------- dirRev
template <typename DXQY>
inline int Boundary<DXQY>::dirRev(const int dir) const
/* dirRev returns the reverse direction of non-zero velocity dir.
 *
 * dir : lattice direction
 */
{
    return DXQY::reverseDirection(dir);
}


#endif // LBBOUNDARY_H

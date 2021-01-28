#ifndef LBINLETOUTLET_H
#define LBINLETOUTLET_H

#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBnodes.h"
#include "LBfield.h"


// *****************************************************
//                        INLET
// *****************************************************
template <typename DXQY>
class Inlet
{
public:
    Inlet(const std::vector<int> boundaryNodes) : boundaryNodes_(boundaryNodes) {}
private:
    std::vector<int> boundaryNodes_;
};


template<>
class Inlet<D2Q9>
{
public:
    Inlet(const std::vector<int> boundaryNodes) : boundaryNodes_(boundaryNodes) {}
    template<typename T>
    void apply(const int &fieldNo, LbField<D2Q9> &f, const Grid<D2Q9> &grid, const lbBase_t & rho, const T &F);
private:
    std::vector<int> boundaryNodes_;
};


template<typename T>
void Inlet<D2Q9>::apply(const int &fieldNo, LbField<D2Q9> &f, const Grid<D2Q9> &grid, const lbBase_t & rho, const T &F)
{
    for (auto nodeNo: boundaryNodes_) {
        lbBase_t jx, jy;
        jx = rho -f(fieldNo, 2, nodeNo) - 2*f(fieldNo, 3, nodeNo) - 2*f(fieldNo, 4, nodeNo)
             - 2*f(fieldNo, 5, nodeNo) - f(fieldNo, 6, nodeNo) - f(fieldNo, 8, nodeNo);
        jy = 1.5*( f(fieldNo, 2, nodeNo) - f(fieldNo, 6, nodeNo) );
        // Set unknown distributions
        f(fieldNo, 0, nodeNo) = f(fieldNo, 4, nodeNo) + 2*D2Q9::c2Inv*D2Q9::w[0]*jx;
        f(fieldNo, 1, nodeNo) = f(fieldNo, 5, nodeNo) + 2*D2Q9::c2Inv*D2Q9::w[1]*(jx + jy);
        f(fieldNo, 7, nodeNo) = f(fieldNo, 3, nodeNo) + 2*D2Q9::c2Inv*D2Q9::w[7]*(jx - jy);
    }
}

// *****************************************************
//                        OUTLET
// *****************************************************
template <typename DXQY>
class Outlet
{
public:
    Outlet(const std::vector<int> boundaryNodes) : boundaryNodes_(boundaryNodes) {}
private:
    std::vector<int> boundaryNodes_;
};


template<>
class Outlet<D2Q9>
{
public:
    Outlet(const std::vector<int> boundaryNodes) : boundaryNodes_(boundaryNodes) {}
    template<typename T>
    void apply(const int &fieldNo, LbField<D2Q9> &f, const Grid<D2Q9> &grid, const lbBase_t & rho, const T &F);
private:
    std::vector<int> boundaryNodes_;
};

template <typename T>
void Outlet<D2Q9>::apply(const int &fieldNo, LbField<D2Q9> &f, const Grid<D2Q9> &grid, const lbBase_t & rho, const T &F)
{
    for (auto nodeNo: boundaryNodes_) {
        lbBase_t jx, jy;
        auto fn = f(fieldNo, nodeNo);
        jx =  2*fn[0] + 2*fn[1] + fn[2] + fn[6] + 2*fn[7] + fn[8] - rho;
        jy =  0;

        f(fieldNo, 3, nodeNo) = -fn[0]/3. - fn[1]/3. - 2*fn[2]/3. + fn[6]/3. + 2*fn[7]/3. - fn[8]/6. + rho/6.;
        f(fieldNo, 4, nodeNo) =  fn[0] + 2*D2Q9::c2Inv*D2Q9::w[4]*(-jx);
        f(fieldNo, 5, nodeNo) = -fn[0]/3. + 2*fn[1]/3. + fn[2]/3. - 2*fn[6]/3. - fn[7]/3. - fn[8]/6. + rho/6.;
    }
}




#endif

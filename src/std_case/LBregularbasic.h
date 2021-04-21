#ifndef LBREGULARBASIC_H
#define LBREGULARBASIC_H

#include <vector>
#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBboundary.h"

//  Linear  package
#include "../Eigen/Dense"
#include "../Eigen/SVD"




template <typename DXQY>
class BoundaryRegularizedBasic
{
public:
    BoundaryRegularizedBasic( const ScalarField       & marker, 
                              const ScalarField       & dist, 
                              const VectorField<DXQY> & normals, 
                              const Nodes<DXQY>       & nodes, 
                              const Grid<DXQY>        & grid, 
                              const LBvtk<DXQY>       & vtk );
    void apply();
    
private:
    std::vector<BDCSVD<MatrixXd>> svd_;
    std::vector<int> nKnown_;
    lbBase_t boundaryMass_;
};



template <typename DXQY>
BoundaryRegularizedBasic<DXQY>::BoundaryRegularizedBasic( const ScalarField       & marker, 
                                                          const ScalarField       & dist, 
                                                          const VectorField<DXQY> & normals, 
                                                          const Nodes<DXQY>       & nodes, 
                                                          const Grid<DXQY>        & grid,
                                                          const LBvtk<DXQY>       & vtk  )
{
    // Count number of boundary nodes
    int nBoundary = 0;
    for (int n = vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n)
    {
        if (marker(0, n) > 0)
        {
            nBoundary += 1;
        }
    }
    
    std::cout << "Number of boundary nodes are " << nBoundary << std::endl;
}


#endif
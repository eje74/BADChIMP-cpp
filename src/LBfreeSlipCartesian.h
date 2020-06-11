#ifndef LBFREESLIPCARTESIAN_H
#define LBFREESLIPCARTESIAN_H

#include "LBglobal.h"
#include "LBboundary.h"
#include "LBgrid.h"
#include "LBfield.h"


/*********************************************************
 * class FREESLIPCARTESIAN: class that performes the
 * free slip boundary condition for static walls directed 
 * along a Cartesian direction, and has
 * the Boundary class as a parent.
 *
 * Here we have assumed that boundary solids are included
 * in the data structure so that we can extract the
 * f-values that are propagated from the solid into the
 * wall.
 *
 * FreeSlipCartesian.apply(...) needs to be run straight 
 * after propagation.
 *
 *********************************************************/

template <typename DXQY>
class FreeSlipCartesian : public Boundary<DXQY>
{
public:
    FreeSlipCartesian(const std::vector<int> &normVec, const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid);     
    void apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    
private:    
    const std::vector<int> n_vec;  // Normal vector
    int q_wall;
    std::vector<int> beta_reflection;  // List of reflected beta values
};


template <typename DXQY>
FreeSlipCartesian<DXQY>::FreeSlipCartesian(const std::vector<int> &normVec, const std::vector<int> bndNodes, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
: n_vec(normVec.begin(), normVec.end()),  Boundary<DXQY>(bndNodes, nodes, grid), beta_reflection(DXQY::nQ)
{
    // Set the direction of the free slip wall
    q_wall = this->dirRev(DXQY::c2q(n_vec));

    
    int n_arr[DXQY::nD];
    std::copy(n_vec.begin(),n_vec.end(),n_arr);
    std::cout<<"Boundary normal: ";
    for (int d = 0; d < DXQY::nD; ++d) {
      std::cout<<n_arr[d]<<" ";
    }
    std::cout<<std::endl;
    
    // Setup the reflected direction
    for (int q = 0; q < DXQY::nQ; ++q ) {
      std::vector<int> c_reflection(DXQY::nD,0);
        for (int d = 0; d < DXQY::nD; ++d) {
            c_reflection[d] = DXQY::c(q, d) - 2 * DXQY::cDot(q, n_arr) * n_vec[d];
        }
        beta_reflection[q] = DXQY::c2q(c_reflection);
	if(beta_reflection[q]==-1){
	  std::cout<<"ERROR in FreeSlipCartesian initialization: c2q returns -1 for q = "<<q<<std::endl;
	}
    }
}


template <typename DXQY>
inline void FreeSlipCartesian<DXQY>::apply(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
/* apply : performs the free slip condition, at the bondary nodes.
 *
 * fieldNo : the lB-field number
 * f       : the field object
 * grid    : grid object
 *
 * Use 'this->' to access functions and variables in the parent class, Boundary<DXQY>.
 *
 */
{
    for (int n = 0; n < this->nBoundaryNodes_; ++n) {
        int node = this->nodeNo(n);
        int node_wall = grid.neighbor(q_wall, node);
        
        for (auto beta: this->beta(n))
        {
            f(fieldNo, beta, node) = f(fieldNo, beta_reflection[beta], node_wall);
        }
    }
}


#endif // LBFREESLIPCARTESIAN_H

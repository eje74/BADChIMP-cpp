#ifndef LBDIFFUSION_H
#define LBDIFFUSION_H

#include "../lbsolver/LBglobal.h"
#include "../lbsolver/LBlatticetypes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBboundary.h"

#include <array>
#include<set>


template<int ND>
std::valarray<lbBase_t> dotLowTriVec(const std::valarray<lbBase_t> &mat, const std::valarray<lbBase_t> &vec)
/*
 *  u_i = mat_ij vec_j
 */
{
    std::cout << "No defiend matrix vector product  for " << ND << " dimensions\n"; 
    return std::valarray<lbBase_t>();
}

template<>
std::valarray<lbBase_t> dotLowTriVec<2>(const std::valarray<lbBase_t> &mat, const std::valarray<lbBase_t> &vec)
/*
 *  u_i = mat_ij vec_j
 */
{
    const lbBase_t v0 = mat[0]*vec[0] + mat[1]*vec[1];
    const lbBase_t v1 = mat[1]*vec[0] + mat[2]*vec[1];

    return std::valarray<lbBase_t> {v0, v1};
}

template<>
std::valarray<lbBase_t> dotLowTriVec<3>(const std::valarray<lbBase_t> &mat, const std::valarray<lbBase_t> &vec)
/*
 *  u_i = mat_ij vec_j
 */
{
    const lbBase_t v0 = mat[0]*vec[0] + mat[1]*vec[1] + mat[3]*vec[2];
    const lbBase_t v1 = mat[1]*vec[0] + mat[2]*vec[1] + mat[4]*vec[2];
    const lbBase_t v2 = mat[3]*vec[0] + mat[4]*vec[1] + mat[5]*vec[2];

    return std::valarray<lbBase_t> {v0, v1, v2};
}

template<int ND>
lbBase_t dotDotLowTriVec(const std::valarray<lbBase_t> &mat, const std::valarray<lbBase_t> &v1, const std::valarray<lbBase_t> &v2)
/*
 * c = v1_i mat_ij v2_j  
 * 
 *  - Assuming that mat is symmetric and on an lower triangular form
 * 
 */
{
    std::cout << "No defiend matrix vector vector contraction for " << ND << " dimensions\n"; 
    return 0;
}

template<>
lbBase_t dotDotLowTriVec<2>(const std::valarray<lbBase_t> &mat, const std::valarray<lbBase_t> &v1, const std::valarray<lbBase_t> &v2)
/*
 * c = v1_i mat_ij v2_j  
 * 
 *  - Assuming that mat is symmetric and on an lower triangular form
 * 
 */
{
    const auto mat_v1 = dotLowTriVec<2>(mat, v1);
    const lbBase_t v2_mat_v1 = v2[0]*mat_v1[0] + v2[1]*mat_v1[1];
    return v2_mat_v1;
}

template<>
lbBase_t dotDotLowTriVec<3>(const std::valarray<lbBase_t> &mat, const std::valarray<lbBase_t> &v1, const std::valarray<lbBase_t> &v2)
/*
 * c = v1_i mat_ij v2_j  
 * 
 *  - Assuming that mat is symmetric and on an lower triangular form
 * 
 */
{
    const auto mat_v1 = dotLowTriVec<3>(mat, v1);
    const lbBase_t v2_mat_v1 = v2[0]*mat_v1[0] + v2[1]*mat_v1[1] + v2[2]*mat_v1[2];
    return v2_mat_v1;
}

template<typename DXQY>
ScalarField readSignedDistance(const std::string &attr, LBvtk<DXQY> & vtklb, const Nodes<DXQY> & nodes, const Grid<DXQY> & grid)
/*
 * Reads the attribute "attr" from the vtklb-file.
 * - Ensures that the sd >= 0 for the fluid and that sd < 0 for the solid
 *
 */
{
    ScalarField sd(1, grid.size()); 
    vtklb.toAttribute(attr);
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        lbBase_t val = vtklb.template getScalarAttribute<lbBase_t>();
        // positive values (and zero) are fluid nodes
        val = std::abs(val);
        if (nodes.isSolid(nodeNo))
            val = -val;
        sd(0, nodeNo) = val;
    }
    return sd;
}

template<typename DXQY>
int addPressureBoundary(const std::string &attr, LBvtk<DXQY> & vtklb, Nodes<DXQY> & nodes, const Grid<DXQY> & grid)
/*
 *  Adds pressure boundaries to the nodes-object.
 *  Returns the number of pressure boundaries.
 *
 *  Reads the "attr" attribute from the vtklb-file.
 *  - Assums integer entries:
 *       0 : do nothing
 *      -1 : set to solid boundary
 *     > 0 : set to fluid boundary if it is a fluid node.
 *           set the nodes tag equal to the attribute.
 *    
 */
{
    int maxBoundaryIndicatorLocal = 0;
    vtklb.toAttribute(attr);
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        maxBoundaryIndicatorLocal = std::max(maxBoundaryIndicatorLocal, pInd);
        if (pInd == -1) { // Solid boundary
            nodes.addNodeType(1, nodeNo);
	    nodes.addNodeTag(pInd, nodeNo);
        } 
        else{
	  nodes.addNodeTag(0, nodeNo);
	}
        
    }

    vtklb.toAttribute(attr);
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        
        if (pInd > 0) { // Fluid boundary
	  if (nodes.isFluid(nodeNo)) {
	    for(int q = 0; q< DXQY::nQ; q++){
	      if(nodes.getTag(grid.neighbor(q, nodeNo)) == -1){
		nodes.addNodeType(2, nodeNo);
		nodes.addNodeTag(pInd, nodeNo);
		break;
	      }
	    }
	  }
        }
        
    }
    
    // Maximum boundary indicator 
    int maxBoundaryIndicator;
    MPI_Allreduce(&maxBoundaryIndicatorLocal, &maxBoundaryIndicator, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);  
    return maxBoundaryIndicator;  
}



template<typename DXQY, typename T>
VectorField<DXQY> calcLaplaceForcing(const int fieldNo, const LbField<DXQY> & f, const lbBase_t tau, const T & bulkNodes) 
{
    VectorField<DXQY> ret(1, f.getNumNodes());
    const lbBase_t tauInv = 1.0/tau;

    for (auto & nodeNo: bulkNodes) {
        ret.set(0, nodeNo) = (DXQY::c2Inv * tauInv) * DXQY::qSumC(f(fieldNo, nodeNo));
    }    
    return ret;
}


template<typename DXQY>
class laplaceBoundary
{
public:
    laplaceBoundary() {}

    template<typename T>
    laplaceBoundary(const T & fluidNodes, const ScalarField &sd, Nodes<DXQY> nodes, const Grid<DXQY> & grid) {
        const auto boundaryNodesTmp = findBoundaryNodes(fluidNodes, sd, nodes, grid);
        boundary_ = Boundary<DXQY>(boundaryNodesTmp[0], nodes, grid);

        int cnt = 0;
        for (const auto & b : boundary_()) {
            const int nodeNo = b.nodeNo();
            const auto normVec = normalFromSignedDistance(nodeNo, sd, grid);
            norms_.push_back(normVec);
            const auto tmp = findWallDirPos(b, sd, normVec, nodes, grid);
            distWall_.push_back(tmp.s);
            dirWallNodes_.push_back(tmp.beta);
            bulkNodes_.push_back(boundaryNodesTmp[1][cnt]);
            boundaryType_.push_back(nodes.getTag(nodeNo));
            const int bulkNode = boundaryNodesTmp[1][cnt];
            std::valarray<lbBase_t> rb(DXQY::nD);
            for (int i = 0; i < DXQY::nD; ++i)
                rb[i] = grid.pos(bulkNode, i) - grid.pos(nodeNo, i);
            rBulk_.push_back(rb); 

            cnt++;
        }
    }

    const Boundary<DXQY> & getBoundary() const {return boundary_;}
    const std::vector<std::valarray<lbBase_t>> & getNorms() const {return norms_;}
    //const std::vector<lbBase_t> & getS() const {return s;}
    const std::vector<int> & getDirWallNodes() const {return dirWallNodes_;}
    /*const std::vector<int> &getBulkNodes() const {return bulkNodes_;}
    const std::vector<int> &getBoundaryType() const {return boundaryType_;} */


    ScalarField readSignedDistance(const std::string & attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> & nodes, const Grid<DXQY> & grid)
    {
        vtklb.toAttribute(attributeName);
        ScalarField sd(1, grid.size());
        for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
        {
            lbBase_t val = vtklb.template getScalarAttribute<lbBase_t>();
            // positive values (and zero) are fluid nodes
            val = std::abs(val);
            if (nodes.isSolid(nodeNo))
                val = -val;
            sd(0, nodeNo) = val;
        }

        return sd;
    }

    std::vector<int>  findBulkFluidNeighbor(const int nodeNo, const ScalarField &sd, const Nodes<DXQY> &nodes,  const Grid<DXQY> &grid)
    /* Search the neighborhood of node 'nodeNo' for fluid neighbors. It will try to 
    * find neighbors as close to the surface normal direction as possilbe 
    *
    * input  : node number, signed distance field, nodes-object, grid-object
    * output : {bulkfluid neighbor, fluid boundary neighbor} -1:not found.
    */
    {
        std::set<int> neigbors;
        for (int q=0; q < DXQY::nQ-1; ++q)
            if (nodes.isFluid(grid.neighbor(q, nodeNo)) && nodes.isMyRank(grid.neighbor(q, nodeNo)))
                neigbors.insert(q);

        int bulkFluid = -1;
        lbBase_t sdBulk = -1.0;
        int boundaryFluid = -1;
        lbBase_t sdBoundary = -1.0;

        for (const auto &q : neigbors) {
            const int neigNo = grid.neighbor(q, nodeNo);
            const lbBase_t sdTest = sd(0, neigNo) / DXQY::cNorm[q];
            if (nodes.isBulkFluid(neigNo)) {
                if (sdTest > sdBulk) {
                    bulkFluid = neigNo;
                    sdBulk = sdTest;
                }
            } else if (nodes.isFluidBoundary(nodeNo)) {
                if (sdTest > sdBoundary) {
                    boundaryFluid = neigNo;
                    sdBoundary = sdTest;
                }
            }
        }    

        return {bulkFluid, boundaryFluid};
    }

    template<typename T>
    std::vector<std::vector<int>> findBoundaryNodes(const T & fluidNodes, const ScalarField &sd, Nodes<DXQY> nodes, const Grid<DXQY> & grid)
    /* Returns a list of bounadry nodes with a list of accompanying bulk node. 
    * The boundary node list is sorted so that a node has an bulk node give that 
    * we can use already calculated boundary nodes as bulk nodes.
    * 
    * input  : list of all fluid nodes in the system, a signed distance function, Nodes-object and Grid-object
    * output : {fluid boundary list, fluid boundaries bulk fluid list}
    */
    {
        std::vector<int> fluidBoundaryNodesRest;
        for (const auto &nodeNo : fluidNodes) 
            if (nodes.isFluidBoundary(nodeNo))
                fluidBoundaryNodesRest.push_back(nodeNo);

        std::vector<int> fluidBoundaryNodes;
        std::vector<int> fluidBoundaryBulkNeighbors;

        while (!fluidBoundaryNodesRest.empty()) {
            std::vector<int> fluidBoundaryNodesRestTmp;
            std::vector<int> setToBulkFluid;
            for (const auto & nodeNo: fluidBoundaryNodesRest) {
                const auto neigNodes = findBulkFluidNeighbor(nodeNo, sd, nodes, grid);
                const int bulkNeig = neigNodes[0];
                const int boundaryNeig = neigNodes[1];
                if (bulkNeig != -1) {
                    fluidBoundaryNodes.push_back(nodeNo);
                    fluidBoundaryBulkNeighbors.push_back(bulkNeig);
                    setToBulkFluid.push_back(nodeNo);
                } 
                else if (boundaryNeig != -1) {
                    fluidBoundaryNodesRestTmp.push_back(nodeNo);
                }
                else {
                    std::cout << "ERROR when trying to find neighboring nodes\n";
                    std::cout << "  It my be that you have a small isolated region in you system\n";
                    std::cout << "  where it is impossible to defined bulk fluid neighbors.\n";
                    exit(1);
                }
            }
            for (const auto & nodeNo: setToBulkFluid) {
                nodes.addNodeType(3, nodeNo);
            }
            fluidBoundaryNodesRest.swap(fluidBoundaryNodesRestTmp);
        }
        return {fluidBoundaryNodes, fluidBoundaryBulkNeighbors};
    }

    std::valarray<lbBase_t> normalFromSignedDistance(const int nodeNo, const ScalarField &sd, const Grid<DXQY> &grid)
    {
        // Setup the sd values of neighborhood of nodeNo
        std::vector<lbBase_t> sdNeigs(DXQY::nQ);
        int cnt = 0; 
        for (const auto & n: grid.neighbor(nodeNo)) {
            sdNeigs[cnt] = sd(0, n);
            cnt++;
        }
        // Gradient and normal
        std::valarray<lbBase_t> grad = DXQY::grad(sdNeigs);
        auto norm = std::sqrt(DXQY::dot(grad, grad));
        if (norm < lbBaseEps) 
            norm = 1.0;
        grad /= norm;

        return grad;
    }

    template<typename T>
    const auto findWallDirPos(const BoundaryNode<DXQY> & bndNode, const ScalarField &sd, const T & normVec, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid)
    {
        class wallPosRet
        {
        public:
            wallPosRet(): beta(-1), nodeNo(-1), s(-1) {}
            int beta;
            int nodeNo;
            lbBase_t s; 

        } ret;
        lbBase_t c_dot_norm = -100;
        int beta = -1;
        for (const auto &q: bndNode.unknown()) {
            const lbBase_t c_dot_norm_try = DXQY::cDotRef(q, normVec) / DXQY::cNorm[q];
            if (c_dot_norm_try > c_dot_norm) {
                c_dot_norm = c_dot_norm_try;
                beta = DXQY::reverseDirection(q);
            }
        }
	
        if (beta < 0) {
	  std::cout << "ERROR in findWallDirPos(): Could not resolve the wall dir. Unknown size: " << bndNode.unknown().size()
		    << "[Type, Tag] = ["<< nodes.getType(bndNode.nodeNo()) << "," << nodes.getTag(bndNode.nodeNo()) <<"]"<< std::endl; 
            //exit(1);
        }

        ret.beta = beta;
        ret.nodeNo = bndNode.nodeNo();
        const int nodeNo = bndNode.nodeNo();
        const int neigNo = grid.neighbor(beta, nodeNo);
        ret.s = sd(0, nodeNo) / (sd(0, nodeNo) - sd(0, neigNo));

        return ret;
    }

    template<typename T1, typename T2, typename T3>
    void addPressureNodeData(const T1 & pressureNodes, const T2 & dists, const T3 & normals)
    {
        for (int bulkNode = 0; bulkNode < boundary_.size(); ++bulkNode) {
            const int nodeNo = boundary_.nodeNo(bulkNode);
            if (boundaryType_[bulkNode] > 0) {
                for (size_t n=0; n < pressureNodes.size(); ++n) { // Find the node in the input-data
                    if (pressureNodes[n] == nodeNo) {
                        distWall_[bulkNode] = dists[n];
                        norms_[bulkNode] = normals[n];
                    }
                }
            }
        }
    }

    template<typename T>
    lbBase_t norm(const T &v) const 
    /*
     * |v| = sqrt(v_i v_i)
     */
    {
        const lbBase_t v2 = DXQY::dot(v, v);
        return std::sqrt(v2);
    }

    template<typename T>
    auto calcMacroValues(const T &f) 
    /*
     *   ret.phi   = <lbBase_t> sum_a f_a
     *   ret.j_i   = <valarray> sum_a c_ai f_a
     *   ret_pi_ij = <valarray> sum_a (c_aic_aj - cs^2 delta_{ij}) f_a
     * 
     *  - pi is returned on a lower triangular form 
     */
    {
        struct ret
        {
        public:
            ret(const std::valarray<lbBase_t> &piInput, const std::valarray<lbBase_t> &jInput, const lbBase_t phiInput) :
            pi(piInput), j(jInput), phi(phiInput) 
            {}
            std::valarray<lbBase_t> pi;
            std::valarray<lbBase_t> j;
            lbBase_t phi;
        };  
        const lbBase_t phi = DXQY::qSum(f);
        const std::valarray<lbBase_t> j=DXQY::qSumC(f);
        const std::valarray<lbBase_t> pi = DXQY::qSumCCLowTri(f) - phi*DXQY::c2*DXQY::deltaLowTri();

        return ret(pi, j, phi);
    }

    lbBase_t caldDist(const int alpha, const lbBase_t phi, const std::valarray<lbBase_t> & j, const std::valarray<lbBase_t> Pi) const
    /*
     *   Calculates the lb-distribution in the alpha direction using the macroscopic values:
     *
     *      f_a = w_a*(rho + c_ai j_i/cs^2 + Q_aij Pi_ij/2cs^4)
     * 
     *   - Assuming that Pi is on a lower triangular form.
     */
    {
        const lbBase_t w = DXQY::w[alpha];
        const auto c = DXQY::cValarray(alpha);

        return w*(phi + DXQY::c2Inv*DXQY::dot(c, j) + DXQY::c4Inv0_5*(dotDotLowTriVec<DXQY::nD>(Pi, c, c) - DXQY::c2*DXQY::traceLowTri(Pi)));
    }

    std::valarray<lbBase_t> jVecExtrapolate(const lbBase_t tau, const std::valarray<lbBase_t> &r_bulk, const std::valarray<lbBase_t> &j_bulk, const std::valarray<lbBase_t> & Pi_bulk) const
    /*
     *   Approximation of the j(0)_i from node values at position r
     *   
     *   j(0)_i = j_i + Pi_ijr_j/(2tau-1)cs^2
     * 
     *   - r, j_i and Pi are r_bulk, j_bulk and Pi_bulk, respectively.
     *  
     */
    {
        std::valarray<lbBase_t> j_bnd = j_bulk + dotLowTriVec<DXQY::nD>(Pi_bulk, r_bulk)/((2*tau - 1)*DXQY::c2);
        return j_bnd;
    }

    std::valarray<lbBase_t> jVecInterpolate(const lbBase_t tau, const std::valarray<lbBase_t> &rWall,const std::valarray<lbBase_t> &rBulk, const std::valarray<lbBase_t> &jWall, const std::valarray<lbBase_t> &jBulk, const std::valarray<lbBase_t> & PiBulk) const
    /*
     *  Apprixmation of j(0) using two node locations w_i (rWall) and  b_i (rBulk)
     *
     *   j(0)_i = (bj_{wi} + wj_{bi})/(b+w) + (bw_j+wb_j)Pi_{ij}/(2tau-1)(b + w)c2^2
     * 
     *  - Where w_i, b_i, j_w, j_b and Pi are rWall, rBulk, jWall, jBulk and PiBulk, respectively
     *  - b = sqrt(b_ib_i)
     *  - w = sqrt(w_iw_i)
     *
     */
    {
        const lbBase_t rwNorm = norm(rWall);
        const lbBase_t rbNorm = norm(rBulk);

        const std::valarray<lbBase_t> jFirstOrder = rbNorm*jWall + rwNorm*jBulk;
        const std::valarray<lbBase_t> jSecondOrder = dotLowTriVec<DXQY::nD>(PiBulk, rbNorm*rWall + rwNorm*rBulk)/((2*tau - 1)*DXQY::c2);
        return (jFirstOrder + jSecondOrder) / (rwNorm + rbNorm);
    }

    std::valarray<lbBase_t> jVecWall(const lbBase_t tau, const std::valarray<lbBase_t> &wallNorm, const std::valarray<lbBase_t> &rWall,const std::valarray<lbBase_t> &rBulk, const std::valarray<lbBase_t> &jBulk, const std::valarray<lbBase_t> & PiBulk) const
    /*
     *   Approxiamtion of j(0) assuming zero flux at the wall. grad\phi\cdot n = 0
     *
     *    j(0)_i = jb_i + (jsn - jb_i n_i)*n_i
     * 
     *    - jb_i : jVecExtrabolate (from tbulk values) 
     *     -jsn = bn_ij_{wi} + wn_ij_{bi})/(b+w) + (bw_j+wb_j)n_iPi_{ij}/(2tau-1)(b + w)c2^2, with n_ij_{wi} = 0 
     */
    {
        const lbBase_t rwNorm = norm(rWall);
        const lbBase_t rbNorm = norm(rBulk);
        const lbBase_t jBulkDotNorm = DXQY::dot(jBulk, wallNorm);

        // Normal direction
        const lbBase_t jNormFirstOrder = rwNorm*jBulkDotNorm;
        const lbBase_t jNormSecondOrder = dotDotLowTriVec<DXQY::nD>(PiBulk, wallNorm, rbNorm*rWall + rwNorm*rBulk)/((2*tau - 1)*DXQY::c2);
        const lbBase_t jWallNorm = (jNormFirstOrder + jNormSecondOrder) / (rwNorm + rbNorm);
        
        // Tangential direction
        const std::valarray<lbBase_t> jFromBulk = jVecExtrapolate(tau, rBulk, jBulk, PiBulk);
        const lbBase_t jBulkNorm = DXQY::dot(wallNorm, jFromBulk);

        return jFromBulk + (jWallNorm - jBulkNorm) * wallNorm;
    }

    lbBase_t rhoPressureBoundary(const lbBase_t tau, const std::valarray<lbBase_t> &rWall, const std::valarray<lbBase_t> &rBulk, const lbBase_t phiWall, const lbBase_t phiBulk, const std::valarray<lbBase_t> &jNode, const std::valarray<lbBase_t> & PiNode) const
    /*
     * phi(0) = (b phi_w + w phi)/(w + b) + (b w_i + w b_i)j_i/(b+w)tau cs^2 - (bw_iw_j + wb_ib_j)Pi_{ij}/2(b+w)(2tau-1)tau cs^4  
     *
     */
    {
        const lbBase_t rwNorm = norm(rWall);
        const lbBase_t rbNorm = norm(rBulk);

        // Normal direction
        const lbBase_t phiZerothOrder = rbNorm*phiWall + rwNorm*phiBulk;
        const std::valarray<lbBase_t> v1 = rbNorm*rWall + rwNorm*rBulk;
        const lbBase_t phiFirstOrder = DXQY::dot(v1, jNode)/(tau*DXQY::c2);
        const lbBase_t bbPi = dotDotLowTriVec<DXQY::nD>(PiNode, rBulk, rBulk);
        const lbBase_t wwPi = dotDotLowTriVec<DXQY::nD>(PiNode, rWall, rWall);
        const lbBase_t phiSecondOrder = -(rwNorm*bbPi + rbNorm*wwPi)/(2*(2*tau - 1)*tau*DXQY::c4);

        return (phiZerothOrder + phiFirstOrder + phiSecondOrder)/(rwNorm + rbNorm);
    }

    template<typename T>
    lbBase_t rhoPressureSolidWall(const T &f, const std::valarray<lbBase_t> &j, const std::valarray<lbBase_t> &Pi, const BoundaryNode<DXQY> bndNode) const
    /*
     *  phi(0) = f0 + \sum_gamma(fgamma + fgamma_rev) + sum_beta(2*fbeta_rev + 2w_beta c_beta_i j_i/cs^2) + sum_delta(w_delta Q_delta_{ij}Pi_{ij}/cs^4) / (1-sum_delta 2w_delta)
     *
     */
    {
        lbBase_t gSum = f[DXQY::nQNonZero_];
        for (const auto & val: bndNode.gamma())
            gSum += f[val];
        for (const auto & val: bndNode.gammaRev())
            gSum += f[val];
        for (const auto & val: bndNode.betaRev())
            gSum += 2*f[val];

        lbBase_t jSum = 0;
        for (const auto & a: bndNode.beta())
            jSum += DXQY::w[a]*DXQY::cDotRef(a, j);
        jSum *= 2*DXQY::c2Inv; 

        lbBase_t piSum = 0;
        lbBase_t wSum = 0;
        for (const auto & a: bndNode.delta()) {
            const auto c = DXQY::cValarray(a);
            piSum += DXQY::w[a]*(dotDotLowTriVec<DXQY::nD>(Pi, c, c) - DXQY::c2*DXQY::traceLowTri(Pi));
            wSum += DXQY::w[a];
        }
        piSum *= DXQY::c4Inv;
        
        return (gSum + jSum + piSum)/(1 - 2*wSum);
    }

    std::valarray<lbBase_t> pressureBoundary(const lbBase_t phiWall, const lbBase_t tau, const std::valarray<lbBase_t> &rw, const std::valarray<lbBase_t> &rb, const lbBase_t phiBulk, const std::valarray<lbBase_t> &jBulk, const std::valarray<lbBase_t> &piBulk, const BoundaryNode<DXQY> &boundaryNode) const
    {
        const std::valarray<lbBase_t> j = jVecExtrapolate(tau, rb, jBulk, piBulk);
        const lbBase_t phi = rhoPressureBoundary(tau, rw, rb, phiWall, phiBulk, j, piBulk);

        std::valarray<lbBase_t> fUnknown(DXQY::nQ);

        for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {
            fUnknown[alpha] = caldDist(alpha, phi, j, piBulk);
        }

        return fUnknown;
    }

    std::valarray<lbBase_t> solidBoudary(const lbBase_t tau, const std::valarray<lbBase_t> &f, const std::valarray<lbBase_t> &wn, const std::valarray<lbBase_t> &rw, const std::valarray<lbBase_t> &rb, const std::valarray<lbBase_t> &jBulk, const std::valarray<lbBase_t> &piBulk, const BoundaryNode<DXQY> &boundaryNode) const
    {
        const std::valarray<lbBase_t> j = jVecWall(tau, wn, rw, rb, jBulk, piBulk);
        const lbBase_t phi = rhoPressureSolidWall(f, j, piBulk, boundaryNode);

        std::valarray<lbBase_t> fUnknown(DXQY::nQ);

        for (int alpha = 0; alpha < DXQY::nQ; ++alpha) {
            fUnknown[alpha] = caldDist(alpha, phi, j, piBulk);
        }

        return fUnknown;
    }


  void apply(const int fieldNum, LbField<DXQY> & f, const lbBase_t tau, const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, ScalarField & applyBnd)
    {
        int count = 0;
        for (const auto & b: boundary_()) {
            const int nodeNo = b.nodeNo();
            const int bulkNodeNo = bulkNodes_[count];
            const int bt = boundaryType_[count];
            const std::valarray<lbBase_t> fBulkNode = f(fieldNum, bulkNodeNo);
            const std::valarray<lbBase_t> rb = rBulk_[count];
            const auto macVal = calcMacroValues(fBulkNode);
            if (macVal.phi != macVal.phi) {
                //std::cout << bt << std::endl;
                std::cout << bulkNodeNo << std::endl;
                for (const auto &val: fBulkNode) {
                    
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
            if ( bt == 0) { // Solid wall boundary
                const auto wn = norms_[count];
                const std::valarray<lbBase_t> rw = DXQY::cValarray(dirWallNodes_[count])*distWall_[count];
                const std::valarray<lbBase_t> fNode = f(fieldNum, nodeNo);
                const std::valarray<lbBase_t> fNew = solidBoudary(tau, fNode, wn, rw, rb, macVal.j, macVal.pi, b);

                f.set(fieldNum, nodeNo) = fNew;
//                for (const auto & alpha: b.unknown())
//                    f(fieldNum, alpha, nodeNo) = fNew[alpha];

            } else if (bt > 0) { // Pressure boundary   
                const lbBase_t phiWall = bt == fieldNum+1 ? 1 : 0;
                const std::valarray<lbBase_t> rw = -norms_[count]*distWall_[count];

		applyBnd(0, nodeNo) = bt;
		
                const std::valarray<lbBase_t> fNew = pressureBoundary(phiWall, tau, rw, rb, macVal.phi, macVal.j, macVal.pi, b);
		//if(std::abs(DXQY::qSum(fNew) - phiWall) > 1e-10)
		//  std::cout<<"ERROR: mismatch in boundary apply. qSum(fNew) = " << DXQY::qSum(fNew) << std::endl;
		  
		
                f.set(fieldNum, nodeNo) = fNew;
//                for (const auto & alpha: b.unknown())
//                    f(fieldNum, alpha, nodeNo) = fNew[alpha];

		//if(std::abs(DXQY::qSum(f(fieldNum, nodeNo)) - phiWall) > 1e-10)
		//  std::cout<<"ERROR: mismatch in boundary apply. qSum(f) = " << DXQY::qSum(f(fieldNum, nodeNo)) << std::endl;

            } else {
                std::cout << "ERROR in boundary condition. Unrecognized boundary type : " << bt << std::endl;
                exit(0);   
            }

            count++;
        }
    }

private:
    Boundary<DXQY> boundary_;
    std::vector<std::valarray<lbBase_t>> norms_;
    std::vector<std::valarray<lbBase_t>> rBulk_;
    std::vector<lbBase_t> distWall_;
    std::vector<int> dirWallNodes_;
    std::vector<int> bulkNodes_;
    std::vector<int> boundaryType_;
};






//=====================================================================================
//
//                         D I F F U S I O N S O L V E R
//
//=====================================================================================
template<typename DXQY>
class DiffusionSolver
{
public:

    // Bulk diffusion
    DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);

    int maxBoundaryIndicator() const {return maxPressureInidcator_;}

    std::valarray<lbBase_t> setF(const lbBase_t & rho) const;

    template<typename T>
    std::valarray<lbBase_t> omegaBGK(const T & f, const lbBase_t & rho) const;

    lbBase_t diffusionCoefficient() const;

    // Geometry
    std::vector<int> findBulkNodes(LBvtk<DXQY> &vtk, const Nodes<DXQY> &nodes);

    // Boundary condition
    void applyBoundaryCondition(LbField<DXQY> &f, const Grid<DXQY> &grid) const;
    void applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const;

    // Set forcing
    template<typename T>
    VectorField<DXQY> getForcing(const int fieldNo, const LbField<DXQY> & f, const T & bulkNodes) const;

    // collision
    std::valarray<lbBase_t> collision(const std::valarray<lbBase_t> &fNode) const;
    std::valarray<lbBase_t> collision(const std::valarray<lbBase_t> &fNode, lbBase_t &rho) const;

    // Write forces and pressures
    template<typename T>
    void writeFieldsToFile(const LbField<DXQY> & f, const ScalarField &rho, const T & bulkNodes);

    // Help function BEGIN
    Boundary<DXQY> getWallBoundary() {return wallBoundary_.bnd;}
    Boundary<DXQY> getWallPressureBoundaryNodes() {return wallPressureBoundary_.bnd;}
    Boundary<DXQY> getPressureBoundaryNodes() {return pressureBoundary_.bnd;}
    std::vector<std::valarray<lbBase_t>> getWallNormals() {return wallBoundary_.normals;}
    std::vector<int> getWallNeighbors() {return wallBoundary_.neighbors;}
    // Help function END
    std::valarray<lbBase_t> getNormals(const int nodeNo) {return normalVector(0, nodeNo);}
    lbBase_t getSignedDistance(const int nodeNo) {return signedDistance(0, nodeNo);}
//    void fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    ScalarField fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);

    // Boundary diffusion
    ScalarField setupBoundaryNodesTest(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);


private:
    // Boundary diffusion
    void setupBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    template<typename T>
    std::vector<T> readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
    template<typename T>
    std::vector<std::valarray<T>> readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid);
 
 
    // Definition of local structures  
    struct DiffusionBoundaryNodes
    {
        std::vector<std::valarray<lbBase_t>> normals;
        std::vector<lbBase_t> qs;
        std::vector<int> indicators;
        std::vector<int> neighbors;
        Boundary<DXQY> bnd;
    };  
    const int size_;
    const lbBase_t tau_;
    const lbBase_t tauInv_;
    const std::valarray<lbBase_t> w_;
    int maxPressureInidcator_;

    DiffusionBoundaryNodes wallBoundary_;
    DiffusionBoundaryNodes pressureBoundary_;
    DiffusionBoundaryNodes wallPressureBoundary_;

    // Test varaibles
    VectorField<DXQY> normalVector;
    ScalarField signedDistance;
};

//                               DiffusionSolver
//----------------------------------------------------------------------------------- DiffusionSolver
//                               DiffusionSolver
//----------------------------------------------------------------------------------- fillBoundaryNodes
template<typename DXQY>
//void DiffusionSolver<DXQY>::fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
ScalarField DiffusionSolver<DXQY>::fillBoundaryNodes(LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
{

    // All nodes
    auto sd = readScalarValues<lbBase_t>("signed_distance", vtk, nodes, grid);
    for (int nodeNo=0; nodeNo<grid.size(); ++nodeNo)
    {
        signedDistance(0, nodeNo) = sd[nodeNo];
        if ( nodes.isFluid(nodeNo) )
        {
            std::vector<lbBase_t> sdNeigs(DXQY::nQ);
            int  cnt = 0;
            for (const auto &neigNo : grid.neighbor(nodeNo)) {
                sdNeigs[cnt] = sd[neigNo];
                cnt += 1;
            }
            std::valarray<lbBase_t> grad = DXQY::grad(sdNeigs);
            auto norm = std::sqrt(DXQY::dot(grad, grad));
            if (norm < lbBaseEps) 
                norm = 1.0;
            grad /= norm;
            normalVector.set(0, nodeNo) = grad;
        }
    }

    // Pressure and wall-pressure boundary normal
    // Read the boundary normal from the vtklb file
    std::vector<std::string> cartExt{"_x", "_y", "_z"};
    for (int d = 0; d < DXQY::nD; ++d) 
    { 
        vtk.toAttribute("boundary_normal" + cartExt[d]);
        for (int n=0; n < vtk.numSubsetEntries(); ++n) {
            const auto att = vtk.template getSubsetAttribure<lbBase_t>();
            normalVector(0, d, att.nodeNo) = att.val; 
            //std::cout << "node[" << d << "] = " << att.nodeNo << " " << att.val << std::endl;
        }
    }

    /* std::valarray<lbBase_t> givenNorm{1.0, 0.0};
    for (int bndNo = 0; bndNo<pressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = pressureBoundary_.bnd.nodeNo(bndNo);   
        normalVector.set(0, nodeNo) = givenNorm;
    }

    for (int bndNo = 0; bndNo<wallPressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = wallPressureBoundary_.bnd.nodeNo(bndNo);   
        normalVector.set(0, nodeNo) = givenNorm;
    } */


    ScalarField vtkPrint(3, grid.size());

    // Set the bulk neighbors  and qvalues
    std::vector<lbBase_t> qValues(grid.size());
    std::vector<int> bulkNeigh(grid.size());


    // wall nodes
    for (int bndNo = 0; bndNo<wallBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = wallBoundary_.bnd.nodeNo(bndNo);  
        // Find the lattice direction that best matches the normal direction
        auto ncVec = DXQY::cDotAll(normalVector(0, nodeNo));
        lbBase_t maxVal = 0;
        int bulkDir = -1;
        lbBase_t sVal = 0;

        for (int q=0; q < DXQY::nQ-1; ++q)
        {            
            const int neigNo = grid.neighbor(q, nodeNo);
            const int neigRevNo = grid.neighbor(DXQY::reverseDirection(q), nodeNo);
            if ( nodes.isMyRank(neigNo) && nodes.isSolid(neigRevNo) && (!nodes.isDefault(neigRevNo)) && nodes.isFluid(neigNo) && (ncVec[q]/DXQY::cNorm[q] > maxVal) ) {
                bulkDir = q;
                maxVal = ncVec[q]/DXQY::cNorm[q];
                sVal = sd[nodeNo]/(sd[nodeNo] - sd[neigRevNo]);
            }
        }
        if (bulkDir == -1) {
            std::cout << "Could not find a bulk node for the diffusion boundary condition" << std::endl;
	        std::cout << "@ node no: " << nodeNo << ", pos ";
	        for(const auto &i: grid.pos(nodeNo))
	            std::cout << i << " ";
	        std::cout << std::endl;
            vtkPrint(2, nodeNo) = 1;
	        // exit(1);
        } 
        bulkNeigh[nodeNo] = bulkDir;
        qValues[nodeNo] = sVal;

        vtkPrint(0, nodeNo) = qValues[nodeNo];
        vtkPrint(1, nodeNo) = bulkNeigh[nodeNo];
    }


    // pressure nodes
    // lbBase_t s_from_file = 0;
    for (int bndNo = 0; bndNo<pressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = pressureBoundary_.bnd.nodeNo(bndNo);  

        // Find the lattice direction that best matches the normal direction
        auto ncVec = DXQY::cDotAll(normalVector(0, nodeNo));
        lbBase_t maxVal = 0;
        int bulkDir = -1;
        //lbBase_t sVal = 0;
        for (int q=0; q < DXQY::nQ-1; ++q)
        {
            const int neigNo = grid.neighbor(q, nodeNo);
            // const int neigRevNo = grid.neighbor(DXQY::reverseDirection(q), nodeNo);
            if ( nodes.isMyRank(neigNo) && nodes.isFluid(neigNo) && (ncVec[q]/DXQY::cNorm[q] > maxVal) ) {
                bulkDir = q;
                maxVal = ncVec[q]/DXQY::cNorm[q];
                // sVal = s_from_file; // sd[nodeNo]/(sd[nodeNo] - sd[neigRevNo]);
            }
        }
        if (bulkDir == -1) {
            std::cout << "Could not find a bulk node for the diffusion boundary condition" << std::endl;
            vtkPrint(2, nodeNo) = 2;
            // exit(1);
        }
        bulkNeigh[nodeNo] = bulkDir;
        // qValues[nodeNo] = sVal;

        //vtkPrint(0, nodeNo) = qValues[nodeNo];
        vtkPrint(1, nodeNo) = bulkNeigh[nodeNo];
    }


    // wall pressure nodes
    for (int bndNo = 0; bndNo<wallPressureBoundary_.bnd.size(); ++bndNo)
    {
        const int nodeNo = wallPressureBoundary_.bnd.nodeNo(bndNo);  

        // Find the lattice direction that best matches the normal direction
        auto ncVec = DXQY::cDotAll(normalVector(0, nodeNo));
        lbBase_t maxVal = 0;
        int bulkDir = -1;
        // lbBase_t sVal = 0;
        for (int q=0; q < DXQY::nQ-1; ++q)
        {
            const int neigNo = grid.neighbor(q, nodeNo);
            //const int neigRevNo = grid.neighbor(DXQY::reverseDirection(q), nodeNo);
            if ( nodes.isMyRank(neigNo) && nodes.isFluid(neigNo) && (ncVec[q]/DXQY::cNorm[q] > maxVal) ) {
                bulkDir = q;
                maxVal = ncVec[q]/DXQY::cNorm[q];
                // sVal = s_from_file; //  sd[nodeNo]/(sd[nodeNo] - sd[neigRevNo]);
            }
        }
        if (bulkDir == -1) {
            std::cout << "Could not find a bulk node for the diffusion boundary condition" << std::endl;
            vtkPrint(2, nodeNo) = 3;
            // exit(1);
        }
        bulkNeigh[nodeNo] = bulkDir;
        // qValues[nodeNo] = sVal;

        // vtkPrint(0, nodeNo) = qValues[nodeNo];
        vtkPrint(1, nodeNo) = bulkNeigh[nodeNo];
    }


    vtk.toAttribute("boundary_distance");
    for (int n=0; n < vtk.numSubsetEntries(); ++n) {
        const auto att = vtk.template getSubsetAttribure<lbBase_t>();
        qValues[att.nodeNo] = att.val/DXQY::cDotRef(bulkNeigh[att.nodeNo], normalVector(0, att.nodeNo)); 
        vtkPrint(0, att.nodeNo) = qValues[att.nodeNo];
    }

//    auto qvalues = readScalarValues<lbBase_t>("q", vtk, nodes, grid);
    // auto normalvec = readVectorValues<lbBase_t>("normal", vtk, nodes, grid);
//    auto neighvec = readVectorValues<int>("neighbor", vtk, nodes, grid);

    auto fillBoundaryNodes = [this, &qValues, &bulkNeigh](DiffusionBoundaryNodes &diffBnd) 
    {
        for (int n = 0; n < diffBnd.bnd.size(); ++n) {
            const int nodeNo = diffBnd.bnd.nodeNo(n);
            diffBnd.qs.push_back(qValues[nodeNo]);
            diffBnd.normals.push_back(normalVector(0, nodeNo));
//            std::vector<int> nvec(DXQY::nD);
//            for (int d=0; d < DXQY::nD; ++d)  nvec[d] = neighvec[nodeNo][d];
            diffBnd.neighbors.push_back(bulkNeigh[nodeNo]);
        }        
    };

/*    auto fillBoundaryNodes = [this, &qvalues, &neighvec](DiffusionBoundaryNodes &diffBnd) 
    {
        for (int n = 0; n < diffBnd.bnd.size(); ++n) {
            const int nodeNo = diffBnd.bnd.nodeNo(n);
            diffBnd.qs.push_back(qvalues[nodeNo]);
            diffBnd.normals.push_back(normalVector(0, nodeNo));
            std::vector<int> nvec(DXQY::nD);
            for (int d=0; d < DXQY::nD; ++d)  nvec[d] = neighvec[nodeNo][d];
            diffBnd.neighbors.push_back(DXQY::c2q(nvec));
        }        
    }; */

    fillBoundaryNodes(wallBoundary_);
    fillBoundaryNodes(pressureBoundary_);
    fillBoundaryNodes(wallPressureBoundary_);

    return vtkPrint;

}


template<typename DXQY>
DiffusionSolver<DXQY>::DiffusionSolver(const lbBase_t tau, LBvtk<DXQY> & vtk, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
:size_(grid.size()), tau_(tau), tauInv_(1.0/tau), w_(DXQY::w, DXQY::nQ), normalVector(1, grid.size()), signedDistance(1, grid.size()) 
{
    
    setupBoundaryNodes(vtk, nodes, grid);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- setF
template<typename DXQY>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::setF(const lbBase_t & rho) const
//-----------------------------------------------------------------------------------
{
    return w_*rho;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- omegaBGK
template<typename DXQY>
template<typename T>
std::valarray<lbBase_t> DiffusionSolver<DXQY>::omegaBGK(const T & f, const lbBase_t & rho) const
//-----------------------------------------------------------------------------------
{
    return tauInv_*(w_*rho - f);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- diffusionCoefficient
template<typename DXQY>
lbBase_t DiffusionSolver<DXQY>::diffusionCoefficient() const
//-----------------------------------------------------------------------------------
{
    return DXQY::c2 * (tau_ - 0.5);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- findBulkNodes
template<typename DXQY>
std::vector<int> DiffusionSolver<DXQY>::findBulkNodes(LBvtk<DXQY> &vtk, const Nodes<DXQY> &nodes)
{
    std::vector<int> bulkNodes;

    vtk.toAttribute("pressure_boundary");
    for (int n = vtk.beginNodeNo(); n < vtk.endNodeNo(); ++n) 
    {
        const auto pInd = vtk.template getScalarAttribute<int>();
        if ( nodes.isFluid(n) && nodes.isMyRank(n) && (pInd >= 0) )
            bulkNodes.push_back(n);
    }
    return bulkNodes;    
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- applyBoundaryCondition
template<typename DXQY>
void DiffusionSolver<DXQY>::applyBoundaryCondition(LbField<DXQY> &f, const Grid<DXQY> &grid) const
//-----------------------------------------------------------------------------------
{
    for (int n=0; n < f.num_fields(); ++n) {
        applyBoundaryCondition(n, f, grid);
    }
}
//                               DiffusionSolver
//----------------------------------------------------------------------------------- applyBoundaryCondition
template<typename DXQY>
void DiffusionSolver<DXQY>::applyBoundaryCondition(const int fieldNo, LbField<DXQY> &f, const Grid<DXQY> &grid) const
//-----------------------------------------------------------------------------------
{
    auto calcpi = [](const std::valarray<lbBase_t> &fNeig) {
        return DXQY::qSumCCLowTri(fNeig) - DXQY::c2 * DXQY::qSum(fNeig) * DXQY::deltaLowTri();
    };
   auto calcjBukl = [&](const std::valarray<lbBase_t> &fNeig, const int alpha, const std::valarray<lbBase_t> &pi) {
        std::valarray<lbBase_t> ret = DXQY::qSumC(fNeig);
        ret +=  DXQY::contractionLowTriVec(pi, DXQY::c(alpha)) / ( (2*tau_ - 1)*DXQY::c2 );
        return ret;
    };
    auto calcjWall = [&](const std::valarray<lbBase_t> &fNeig, const int alpha, const lbBase_t q, const std::valarray<lbBase_t> &nVec, const std::valarray<lbBase_t> &pi)  {
        const auto jNeig = DXQY::qSumC(fNeig); // \vec{j}(\vec{c}_\beta) in article
        std::valarray<lbBase_t> ret = jNeig;
        ret -= ( DXQY::dot(nVec, jNeig)/(1+q) )*nVec;
        const auto piDotC = DXQY::contractionLowTriVec(pi, DXQY::c(alpha));
        ret += ( piDotC - nVec*DXQY::dot(nVec, piDotC) ) / ( (2*tau_ - 1)*DXQY::c2 );
        return ret;
    };
   auto calcPhi = [&](const std::valarray<lbBase_t> &fNode, const std::vector<int> &gamma, const std::vector<int> &beta, const std::vector<int> &delta, const std::valarray<lbBase_t> &jVec, const std::valarray<lbBase_t> &pi) {
        lbBase_t rhs = 0.0;
        rhs += fNode[DXQY::nQNonZero_];
        for (auto & alpha : gamma) {
            rhs += fNode[alpha];
            rhs += fNode[DXQY::reverseDirection(alpha)];
        }
        for (auto alpha: beta) {
            const lbBase_t cj = DXQY::dot(DXQY::c(alpha), jVec);
            rhs += 2*fNode[DXQY::reverseDirection(alpha)] + 2*w_[alpha]*DXQY::c2Inv*cj;
        }
        lbBase_t lhs = 1.0;
        for (auto alpha: delta) {
            const lbBase_t cPic = DXQY::dot(DXQY::c(alpha),DXQY::contractionLowTriVec(pi, DXQY::c(alpha)));
            const lbBase_t piTrace = DXQY::traceLowTri(pi);
            rhs += w_[alpha]*DXQY::c4Inv*(cPic - DXQY::c2*piTrace);
            lhs -= 2*w_[alpha];
        }
        return rhs/lhs;
    };
    auto calcf = [&](const lbBase_t alpha, const lbBase_t phi, const std::valarray<lbBase_t> &jVec, const std::valarray<lbBase_t> &piMat) {
        const lbBase_t cj = DXQY::dot(DXQY::c(alpha), jVec);
        const lbBase_t ccPi = DXQY::dot(DXQY::c(alpha),DXQY::contractionLowTriVec(piMat, DXQY::c(alpha)));
        const lbBase_t iiPi = DXQY::traceLowTri(piMat);

        return w_[alpha]*(phi + DXQY::c2Inv*cj + DXQY::c4Inv0_5*(ccPi - DXQY::c2*iiPi));
    };  

    // ------------------------------------------------------------------------------ WALL BOUNDARY
    const DiffusionBoundaryNodes &diffBnd = wallBoundary_;
    const Boundary<DXQY> &bnd = diffBnd.bnd;
    for (int bndNo = 0; bndNo < bnd.size(); ++bndNo) 
    {
        const int nodeNo = bnd.nodeNo(bndNo);
        const int alphaNeig = diffBnd.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));
        
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig); 

        // Calculate the first moment
        const auto qNode = diffBnd.qs[bndNo];  // q in article
        const auto nNode = diffBnd.normals[bndNo]; // \vec{n} in article

        auto jNode = calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        auto phiNode = calcPhi(f(fieldNo, nodeNo), bnd.gamma(bndNo), bnd.beta(bndNo), bnd.delta(bndNo), jNode, piNode);

        //
        for (auto &beta : bnd.beta(bndNo)) {
            auto betaHat = bnd.dirRev(beta);
            f(fieldNo, beta, nodeNo) = f(fieldNo, betaHat, nodeNo) + 2*DXQY::c2Inv*w_[beta]*DXQY::dot(DXQY::c(beta), jNode);
        }

        for (auto &delta : bnd.delta(bndNo)) {
            f(fieldNo, delta, nodeNo) = calcf(delta, phiNode, jNode, piNode);
            auto deltaHat = bnd.dirRev(delta);
            f(fieldNo, deltaHat, nodeNo) = calcf(deltaHat, phiNode, jNode, piNode);
        }
    }

    // ------------------------------------------------------------------------------ WALL PRESSURE BOUNDARY 
    const DiffusionBoundaryNodes &diffBndWP = wallPressureBoundary_;
    const Boundary<DXQY> &bndWP = diffBndWP.bnd;
    for (int bndNo = 0; bndNo < bndWP.size(); ++bndNo)
    {
        const int nodeNo = bndWP.nodeNo(bndNo);
        const int alphaNeig = diffBndWP.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig);
        // Calculate the first moment
        const auto qNode = diffBndWP.qs[bndNo]; // q in article
        // const auto nNode = diffBndWP.normals[bndNo]; // \vec{n} in article
        auto jNode = calcjBukl(fNeig, alphaNeig, piNode); // calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        // Calculate the zeroth moment
        const lbBase_t phiWall = 1.0 * (diffBndWP.indicators[bndNo] == (fieldNo + 1));
        const lbBase_t phiBulk = DXQY::qSum(fNeig);
        const lbBase_t ccPi = DXQY::dot(DXQY::c(alphaNeig), DXQY::contractionLowTriVec(piNode, DXQY::c(alphaNeig)));

        const lbBase_t phiNode = (phiWall + qNode * phiBulk) / (1 + qNode) - 0.5 * qNode * DXQY::c4Inv * ccPi / ((2 * tau_ - 1.0) * tau_);

        auto alphaZero = DXQY::nQNonZero_;
        f(fieldNo, alphaZero, nodeNo) = calcf(alphaZero, phiNode, jNode, piNode);
        auto setfsRev = [&](const std::vector<int> &alphaList)
        {
            for (auto &alpha : alphaList)
            {
                auto alphaHat = bndWP.dirRev(alpha);
                f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
            }
        };
        auto setfs = [&](const std::vector<int> &alphaList)
        {
            for (auto &alpha : alphaList)
            {
                f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
            }
        };
        // setfs(bndWP.gamma(bndNo));
        setfs(bndWP.beta(bndNo));
        setfs(bndWP.delta(bndNo));
        setfsRev(bndWP.delta(bndNo));
    }
    
    // ------------------------------------------------------------------------------ PRESSURE BOUNDARY
    const DiffusionBoundaryNodes &diffBndP = pressureBoundary_;
    const Boundary<DXQY>  &bndP = diffBndP.bnd;
    for (int bndNo = 0; bndNo < bndP.size(); ++bndNo)
    {
        const int nodeNo = bndP.nodeNo(bndNo);
        const int alphaNeig = diffBndP.neighbors[bndNo];
        const auto fNeig = f(fieldNo, grid.neighbor(alphaNeig, nodeNo));
        // Calculate the second moment
        const std::valarray<lbBase_t> piNode = calcpi(fNeig);
        // Calculate the first moment
        const auto qNode = diffBndP.qs[bndNo];  // q in article
        //const auto nNode = diffBndP.normals[bndNo]; // \vec{n} in article
        auto jNode = calcjBukl(fNeig, alphaNeig, piNode);//calcjWall(fNeig, alphaNeig, qNode, nNode, piNode);
        // Calculate the zeroth moment
        const lbBase_t phiWall = 1.0 * (diffBndP.indicators[bndNo] == (fieldNo+1));
        const lbBase_t phiBulk = DXQY::qSum(fNeig);
        const lbBase_t ccPi = DXQY::dot(DXQY::c(alphaNeig),DXQY::contractionLowTriVec(piNode, DXQY::c(alphaNeig)));

        const lbBase_t phiNode = (phiWall + qNode*phiBulk)/(1 + qNode) - 0.5*qNode*DXQY::c4Inv*ccPi/((2*tau_-1.0)*tau_);

        auto alphaZero = DXQY::nQNonZero_;
        f(fieldNo, alphaZero, nodeNo) = calcf(alphaZero, phiNode, jNode, piNode);

    /*        auto setfs = [&](const std::vector<int> &alphaList) {
                for (auto & alpha: alphaList) {
                    f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
                    auto alphaHat = bndP.dirRev(alpha);
                    f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
                }
            };
            setfs(bndP.gamma(bndNo));
            setfs(bndP.beta(bndNo));
            setfs(bndP.delta(bndNo));
    */
        auto setfsRev = [&](const std::vector<int> &alphaList) {
                for (auto & alpha: alphaList) {
                    auto alphaHat = bndP.dirRev(alpha);
                    f(fieldNo, alphaHat, nodeNo) = calcf(alphaHat, phiNode, jNode, piNode);
                }
            };
            auto setfs = [&](const std::vector<int> &alphaList) {
                for (auto & alpha: alphaList) {
                    f(fieldNo, alpha, nodeNo) = calcf(alpha, phiNode, jNode, piNode);
                }
            };
            setfs(bndP.beta(bndNo));
            setfs(bndP.delta(bndNo));
            setfsRev(bndP.delta(bndNo));

    } 
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- setupBoundaryNodesTest
template<typename DXQY>
ScalarField DiffusionSolver<DXQY>::setupBoundaryNodesTest(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    std::vector<int> wallBoundaryNodes;
    std::vector<int> pressureBoundaryNodes;
    std::vector<int> wallPressureBoundaryNodes;

    ScalarField ret(1, grid.size());

    // Read pressure_boundary
    int maxPressureInidcator = 0;
    vtklb.toAttribute("pressure_boundary");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        maxPressureInidcator = std::max(maxPressureInidcator, pInd);

        ret(0, nodeNo) = 0;

        if ( (nodes.isFluidBoundary(nodeNo)  || (pInd > 0)) && nodes.isFluid(nodeNo) && nodes.isMyRank(nodeNo) && (pInd >= 0)) {
            auto hasSolidNeighbors = [&nodeNo, &nodes, &grid]() -> bool {
                for (auto neighNo: grid.neighbor(nodeNo)) 
                    if ( nodes.isBulkSolid(neighNo) || nodes.isSolidBoundary(neighNo) ) 
                        return true;
                return false;
            };

            if (hasSolidNeighbors()) {
                if (pInd == 0) {
                    wallBoundaryNodes.push_back(nodeNo);
                } else {
                    wallPressureBoundaryNodes.push_back(nodeNo);
                    wallPressureBoundary_.indicators.push_back(pInd);
                }
            } else {
                if (pInd > 0) {
                    pressureBoundaryNodes.push_back(nodeNo);
                    pressureBoundary_.indicators.push_back(pInd);
                } else {
                    std::cout << "Node " << nodeNo << " on processor rank " << nodes.getRank(nodeNo);
                    std::cout << " is not a wall node and is not assigned a pressure!" << std::endl;
                    //exit(1);
                    ret(0, nodeNo) = 1;
                }
            }
        }  
        /* // Setup the different boundary nodes.
        wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
        pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
        wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid); */  
    }

    return ret;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- setupBoundaryNodes
template<typename DXQY>
void DiffusionSolver<DXQY>::setupBoundaryNodes(LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    std::vector<int> wallBoundaryNodes;
    std::vector<int> pressureBoundaryNodes;
    std::vector<int> wallPressureBoundaryNodes;

    // Read pressure_boundary
    int maxPressureInidcator = 0;
    vtklb.toAttribute("pressure_boundary");
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        const auto pInd = vtklb.template getScalarAttribute<int>();
        maxPressureInidcator = std::max(maxPressureInidcator, pInd);
        if ( (nodes.isFluidBoundary(nodeNo)  || (pInd > 0)) && nodes.isFluid(nodeNo) && nodes.isMyRank(nodeNo) && (pInd >= 0)) {
            auto hasSolidNeighbors = [&nodeNo, &nodes, &grid]() -> bool {
                for (auto neighNo: grid.neighbor(nodeNo)) 
                    if ( nodes.isBulkSolid(neighNo) || nodes.isSolidBoundary(neighNo) ) 
                        return true;
                return false;
            };

            if (hasSolidNeighbors()) {
                if (pInd == 0) {
                    wallBoundaryNodes.push_back(nodeNo);
                } else {
                    wallPressureBoundaryNodes.push_back(nodeNo);
                    wallPressureBoundary_.indicators.push_back(pInd);
                }
            } else {
                if (pInd > 0) {
                    pressureBoundaryNodes.push_back(nodeNo);
                    pressureBoundary_.indicators.push_back(pInd);
                } else {
                    std::cout << "Node " << nodeNo << " on processor rank " << nodes.getRank(nodeNo);
                    std::cout << " is not a wall node and is not assigned a pressure!" << std::endl;
                    exit(1);
                }
            }
        }  
        /* // Setup the different boundary nodes.
        wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
        pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
        wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid); */  
    }
    // Setup the different boundary nodes.
    wallBoundary_.bnd = Boundary<DXQY>(wallBoundaryNodes, nodes, grid);
    pressureBoundary_.bnd =  Boundary<DXQY>(pressureBoundaryNodes, nodes, grid);
    wallPressureBoundary_.bnd = Boundary<DXQY>(wallPressureBoundaryNodes, nodes, grid);        
    maxPressureInidcator_ = maxPressureInidcator;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- readScalarValues
template<typename DXQY>
template<typename T>
std::vector<T> DiffusionSolver<DXQY>::readScalarValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    vtklb.toAttribute(attributeName);
    std::vector<T> ret(grid.size());
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
    {
        ret[nodeNo] = vtklb.template getScalarAttribute<T>();
    }

    return ret;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- readVectorValues
template<typename DXQY>
template<typename T>
std::vector<std::valarray<T>> DiffusionSolver<DXQY>::readVectorValues(const std::string attributeName, LBvtk<DXQY> & vtklb, const Nodes<DXQY> &nodes, const Grid<DXQY> & grid)
//-----------------------------------------------------------------------------------
{
    const std::array<const std::string, 3> index{"x", "y", "z"};
    std::vector<std::valarray<T>> ret(grid.size(), std::valarray<T>(DXQY::nD));

    for (int d = 0; d < DXQY::nD; ++d) {
        vtklb.toAttribute(attributeName + "_" + index[d]);
        for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); nodeNo++) 
        {
            ret[nodeNo][d] = vtklb.template getScalarAttribute<T>();
        }
    }
    return ret;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- getForcing
template<typename DXQY>
template<typename T>
VectorField<DXQY> DiffusionSolver<DXQY>::getForcing(const int fieldNo, const LbField<DXQY> & f, const T & bulkNodes) const
//-----------------------------------------------------------------------------------
{
    VectorField<DXQY> ret(1, f.getNumNodes());

    for (auto & nodeNo: bulkNodes) {
        ret.set(0, nodeNo) = (DXQY::c2Inv * tauInv_) * DXQY::qSumC(f(fieldNo, nodeNo));
    }    
    return ret;
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- collision
template<typename DXQY>    
std::valarray<lbBase_t> DiffusionSolver<DXQY>::collision(const std::valarray<lbBase_t> &fNode) const
//-----------------------------------------------------------------------------------
{
    const lbBase_t rhoNode = calcRho<DXQY>(fNode);
    return fNode + omegaBGK(fNode, rhoNode);
}

//                               DiffusionSolver
//----------------------------------------------------------------------------------- collision
template<typename DXQY>    
std::valarray<lbBase_t> DiffusionSolver<DXQY>::collision(const std::valarray<lbBase_t> &fNode, lbBase_t &rho) const
//-----------------------------------------------------------------------------------
{
    const lbBase_t rhoNode = calcRho<DXQY>(fNode);
    rho = rhoNode;
    return fNode + omegaBGK(fNode, rhoNode);
}

//                               DiffusionSolver
//-----------------------------------------------------------------------------------
template<typename DXQY>  
template<typename T>
void DiffusionSolver<DXQY>::writeFieldsToFile(const LbField<DXQY> & f, const ScalarField &rho, const T & bulkNodes)
{

}

#endif

#ifndef LBGLOBALFORCING_H
#define LBGLOBALFORCING_H

#include "LBglobal.h"
#include "LBlatticetypes.h"
#include "LBfield.h"

template <typename DXQY>
inline lbBase_t calcFluxForceCartDir(const int fieldNo, LbField<DXQY> &f, const std::vector<int> bulkNodes, const int cartDir, const lbBase_t fixedFlux, const int numNodesGlobal)
/* calcFluxForce  : sets uniform force resulting in fixed global flux in a single direction 
 *
 * f              : pointer to lb distribution
 * bulkNodes      : pointer to bulkNodes
 * cartDir        : integer [0,1, or 2] giving the Cartesian direction in which the flux is fixed
 * fixedFlux      : Value of the fixed flux
 * numNodesGlobal : Number of Nodes in the system of which the average is taken   

 */
{
    lbBase_t ret;
    lbBase_t meanfc=0.0;
    lbBase_t meanfcGlobal;
    for (auto nodeNo : bulkNodes) {
      auto fTot = f(fieldNo, nodeNo);
      auto sumfc = DXQY::qSumC(fTot);
      meanfc+=sumfc[cartDir];
    }
    MPI_Allreduce(&meanfc, &meanfcGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    meanfcGlobal/=numNodesGlobal; 
    
    return ret=2*(fixedFlux-meanfcGlobal);
}

template <typename DXQY>
inline lbBase_t calcCapNumbForceCartDir(const int fieldNo, LbField<DXQY> &f, ScalarField &rho, const std::vector<int> bulkNodes, const int cartDir, const lbBase_t sigmaCapNumb, const lbBase_t nu0, const lbBase_t nu1, const int numNodesGlobal)
/* calcCapNumbForceCartDir : sets uniform force resulting in fixed Capillary Number in a single direction 
 *
 * f              : pointer to lb distribution
 * rho0           : pointer to the density field of all fluid components
 * bulkNodes      : pointer to bulkNodes
 * cartDir        : integer [0,1, or 2] giving the Cartesian direction in which the flux is fixed
 * fixedFlux      : Value of the fixed flux
 * nu0            : kinematic viscosity of fluid 0
 * nu1            : kinematic viscosity of fluid 1
 * numNodesGlobal : Number of Nodes in the system of which the average is taken   
 */
{
    lbBase_t ret;
    lbBase_t meanPhi0fc=0.0;
    lbBase_t meanPhi1fc=0.0;
    
    lbBase_t meanPhi0fcGlobal;
    lbBase_t meanPhi1fcGlobal;

    lbBase_t meanPhi0=0.0;
    lbBase_t meanPhi1=0.0;
    lbBase_t meanPhi0Global;
    lbBase_t meanPhi1Global;
	
    for (auto nodeNo : bulkNodes) {
      auto fTot = f(0, nodeNo);
      auto sumfc = DXQY::qSumC(fTot);

      const auto rho0Node = rho(0, nodeNo); 
      const auto rho1Node = rho(1, nodeNo);
      const auto rhoTotNode = rho0Node + rho1Node;
      
      const auto phi0Node = rho0Node/rhoTotNode;
      const auto phi1Node = rho1Node/rhoTotNode;
      
	  
      meanPhi0fc+=phi0Node*sumfc[cartDir];
      meanPhi1fc+=phi1Node*sumfc[cartDir];

      meanPhi0+=phi0Node;
      meanPhi1+=phi1Node;
    }
	
	
	MPI_Allreduce(&meanPhi0fc, &meanPhi0fcGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meanPhi1fc, &meanPhi1fcGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meanPhi0, &meanPhi0Global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&meanPhi1, &meanPhi1Global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	

	
	meanPhi0fcGlobal/=numNodesGlobal;
	meanPhi1fcGlobal/=numNodesGlobal;

	meanPhi0Global/=numNodesGlobal;
	meanPhi1Global/=numNodesGlobal;

	
	
	return ret=2*(sigmaCapNumb - (meanPhi0fcGlobal*nu0+meanPhi1fcGlobal*nu1))/(meanPhi0Global*nu0+meanPhi1Global*nu1);
	
}

#endif

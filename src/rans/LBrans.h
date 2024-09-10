#ifndef LBRANS_H
#define LBRANS_H

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../lbsolver/LBglobal.h"
#include "../io/Input.h"
#include "LBBoundaryInterpolation.h"
#include "LBBoundaryRans.h"



//=====================================================================================
//
//                                      R A N S
//
//=====================================================================================
template<typename DXQY>
class Rans
{
public:
    Rans() {};
  Rans(Input &input, int myRank);

  template <typename T1, typename T2, typename T3>
  void apply( 		  
	     const T1 &f,
	     const lbBase_t& rho,
	     const T2 &u,
	     const lbBase_t& u_sq,
	     const std::valarray<lbBase_t> &cu,
	     const T3 &F,
	     const lbBase_t &source,
	     const T1 &g,
	     const T1 &h,
	     const lbBase_t& rhoKInput,
	     const lbBase_t& rhoEInput);

    inline lbBase_t tau() const {return tau_;}
    inline lbBase_t tauK() const {return tauK_;}
    inline lbBase_t tauE() const {return tauE_;}

    inline lbBase_t sourceK() const {return sourceK_;}
    inline lbBase_t sourceE() const {return sourceE_;}

    inline lbBase_t rhoK() const {return rhoK_;}
    inline lbBase_t rhoE() const {return rhoE_;}

  inline lbBase_t gammaDot() const {return gammaDotTilde_;}

  template <typename T1, typename T2, typename T3>
  void zouHeFixedValueLeftBnd(
		  const T1 &bndNodes,
		  T2 &f,
		  T3 &rho,
		  const lbBase_t fixValue/*,
		  const T3 &grid*/);
  template <typename T1, typename T2/*, typename T3*/>
  void zouHeFixedVelocityLeftBnd(
		  const T1 &bndNodes,
		  T2 &f,
		  const lbBase_t fixValue/*,
		  const T3 &grid*/);
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void zouHeOpenRightBnd(
		  const T1 &bndNodes,
		  T2 &f,
		  const T3 &scalarField,
		  //const lbBase_t fixValue,
		  const T4 &force,
		  const T5 &grid);
  template <typename T1, typename T2, typename T3, typename T4>
  void zouHeFixedValueRightBnd(
		  const T1 &bndNodes,
		  T2 &f,
		  const lbBase_t fixValue,
		  const T3 &force,
		  const T4 &grid);
   template <typename T1, typename T2>
   void copyDistBnd(
			const T1 &bndNodes,
			T2 &f,
      int qDir,
      const Grid<DXQY> &grid); 
   lbBase_t wallShear(lbBase_t uIntrp);
   void solidBnd(
      LbField<DXQY> &f,
      LbField<DXQY> &fTmp,
      ScalarField &rho,
      VectorField<DXQY> &vel,
      ScalarField &viscocity,
      LbField<DXQY> &g,
      LbField<DXQY> &gTmp,
      ScalarField &rhoK, 
      LbField<DXQY> &h,
      LbField<DXQY> &hTmp,
      ScalarField &rhoE,
      std::vector<InterpolationElement> &boundarynodes, 
      Nodes<DXQY> &nodes,
      Grid<DXQY> &grid);
   void solidBndVel(
      LbField<DXQY> &f,
      ScalarField &rho,
      VectorField<DXQY> &vel,
      ScalarField &viscocity,
      LbField<DXQY> &g,
      ScalarField &rhoK, 
      LbField<DXQY> &h,
      ScalarField &rhoE,
      std::vector<InterpolationElement> &boundarynodes, 
      Grid<DXQY> &grid);


   

private:
  //                             Constant input parameters
  //------------------------------------------------------------------------------------- Constant input parameters
    const lbBase_t viscosity0_;
    const lbBase_t tau0_;
    const lbBase_t Cmu_;
    const lbBase_t C1epsilon_;
    const lbBase_t C2epsilon_;
    const lbBase_t sigma0kInv_;
    const lbBase_t sigmakInv_;
    const lbBase_t sigma0epsilonInv_;
    const lbBase_t sigmaepsilonInv_;

    const lbBase_t kappa_;
    const lbBase_t E_;
    const lbBase_t yp_;

    const lbBase_t velNoiseAmplitude_;
  //------------------------------------------------------------------------------------- Constant derived variables
    const lbBase_t X1_;
    const lbBase_t Y1_;
    const lbBase_t Z1_;
    const lbBase_t Z2_;

    const lbBase_t tauK0Term_;
    const lbBase_t tauE0Term_;

    lbBase_t uWallMin_;
    lbBase_t uWallMax_;

    const lbBase_t lowC0_;  // Constants used in the law of the wall calculation
    const lbBase_t lowC1_;

  //                                    parameters
  //------------------------------------------------------------------------------------- Parameters to be calculated
    lbBase_t tau_;
    lbBase_t tauK_;
    lbBase_t tauE_;

    lbBase_t sourceK_;
    lbBase_t sourceE_;

    lbBase_t rhoK_;
    lbBase_t rhoE_;

    lbBase_t gammaDotTilde_;
};

//                                        Rans
//------------------------------------------------------------------------------------- Rans
//NOTE: When initializing const in constructor important that the order of parameters is the same as in declaration list (see above) 
template<typename DXQY>
Rans<DXQY>::Rans(Input &input, int myRank)
  :viscosity0_(input["fluid"]["viscosity"]),
   tau0_(DXQY::c2Inv * input["fluid"]["viscosity"] + 0.5), 
   Cmu_(input["RANS"]["k-epsilonCoef"]["C_mu"]),
   C1epsilon_(input["RANS"]["k-epsilonCoef"]["C_1epsilon"]),
   C2epsilon_(input["RANS"]["k-epsilonCoef"]["C_2epsilon"]),
   sigma0kInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_0k"]),
   sigmakInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_k"]),
   sigma0epsilonInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_0epsilon"]),
   sigmaepsilonInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_epsilon"]),
   kappa_(input["RANS"]["wall"]["kappa"]),
   E_(input["RANS"]["wall"]["E"]),
   yp_(input["RANS"]["wall"]["yp"]),
   velNoiseAmplitude_(input["RANS"]["inlet"]["velNoiseAmplitude"]),
   X1_(Cmu_*DXQY::c2Inv),
   Y1_(0.25*DXQY::c2Inv),
   Z1_(0.25*DXQY::c2Inv*C1epsilon_),
   Z2_(C2epsilon_),
   tauK0Term_((tau0_-0.5)*sigma0kInv_),
   tauE0Term_((tau0_-0.5)*sigma0epsilonInv_),
   uWallMin_(10.0*viscosity0_*std::log(E_*10.0)/(kappa_*yp_)),
   uWallMax_(300.0*viscosity0_*std::log(E_*300.0)/(kappa_*yp_)),
   lowC0_(1.0/kappa_),
   lowC1_(E_*yp_/viscosity0_)

/* Class constructor, sets k-epsilon coefficent constant values from input file
 * 
 * Parameters
 * ----------
 * input : input object containing input file.
 * 
 * Returns
 * -------
 * Rans<DXQY> constructor
 *    DXQY is a lattice type
 */   
{
  if (myRank == 0){
  std::cout << "-------------RANS PARAMETERS-----------------" << std::endl;
  std::cout << "tau0 = " << tau0_ << std::endl;
  std::cout << "Cmu = " << Cmu_ << std::endl;
  std::cout << "C1epsilon = " << C1epsilon_ << std::endl;
  std::cout << "C2epsilon = " << C2epsilon_ << std::endl;
  std::cout << "sigma0k = " << 1./sigma0kInv_ << std::setw(30)<< "sigma0kInv = " << sigma0kInv_ << std::endl;
  std::cout << "sigmak = " << 1./sigmakInv_ << std::setw(30)<< "sigmakInv = " << sigmakInv_ << std::endl;
  std::cout << "sigma0epsilon = " << 1./sigma0epsilonInv_ << std::setw(30)<< "sigma0epsilonInv = " << sigma0epsilonInv_ << std::endl;
  std::cout << "sigmaepsilon = " << 1./ sigmaepsilonInv_ << std::setw(30)<< "sigmaepsilonInv = " << sigmaepsilonInv_ << std::endl;
  std::cout << "kappa = " << kappa_ << std::endl;
  std::cout << "E = " << E_ << std::endl;
  std::cout << "yp = " << yp_ << std::endl;
  std::cout << "X1 = " << X1_ << std::endl;
  std::cout << "Y1 = " << Y1_ << std::endl;
  std::cout << "Z1 = " << Z1_ << std::endl;
  std::cout << "Z2 = " << Z2_ << std::endl;
  std::cout << "uWallMin = " << uWallMin_ << std::endl;
  std::cout << "uWallMax = " << uWallMax_ << std::endl;
  std::cout << std::endl;
  std::cout << "velNoiseAmplitude = " << velNoiseAmplitude_ << std::endl;
  std::cout << std::endl;

  }

  std::srand(0);
}

template <typename DXQY>
template <typename T1, typename T2, typename T3>
void Rans<DXQY>::apply( 		  
		       const T1 &f,
		       const lbBase_t& rho,
		       const T2 &u,
		       const lbBase_t& u_sq,
		       const std::valarray<lbBase_t> &cu,
		       const T3 &F,
		       const lbBase_t &source,
		       const T1 &g,
		       const T1 &h,
		       const lbBase_t& rhoKInput,
		       const lbBase_t& rhoEInput)
/* Returns void: computes object values
 * 
 * 
 * 
 * Parameters
 * ----------
 * f : array-like, size = [DXQY::nQ]
 *     lb distribution
 * 
 * rho : float-like
 *     fluid density
 * 
 * u : array-like, size = [DXQY::nD]
 *     fluid velocity
 * 
 * u_sq : float-like
 *     square of the fluid velocity
 * 
 * cu : valarray, size = [DXQY::nQ]
 *     dot product of ``u`` and the basis velocities
 * 
 * F : array-like, size = [DXQY::nD]
 *     body force
 * 
 * source :  float-like
 *     bulk fluid source
 * 
 * Returns
 * -------
 * void
 *     
 */ 
{

  std::valarray<lbBase_t> feq(DXQY::nQ);
  std::valarray<lbBase_t> strain_rate_tilde(0.0, DXQY::nD * DXQY::nD);
  lbBase_t strain_rate_tilde_square = 0;    

  //                                 Strain rate calculation
  //------------------------------------------------------------------------------------- Strain rate calculation
  for (int q = 0; q < DXQY::nQ; ++q){
    feq[q] = DXQY::w[q]*rho*( 1.0 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5 * (cu[q]*cu[q] - DXQY::c2*u_sq) );        

    for (int i=0; i < DXQY::nD; ++i){
      strain_rate_tilde[i + DXQY::nD*i] -=  (f[q] - feq[q])*( - DXQY::c2);
      for (int j = 0; j < DXQY::nD; ++j){
	      strain_rate_tilde[i + DXQY::nD*j] -=  (f[q] - feq[q])*DXQY::c(q, i)*DXQY::c(q, j);
      }            
    }
  }
  
  for (int i=0; i < DXQY::nD; ++i){
    strain_rate_tilde[i + DXQY::nD*i] -= 0.5*DXQY::c2*source;
    for (int j = 0; j < DXQY::nD; ++j){
      strain_rate_tilde[i + DXQY::nD*j] -= 0.5*(u[i]*F[j] + u[j]*F[i]);
      strain_rate_tilde[i + DXQY::nD*j] -= 0.5*u[i]*u[j]*source;
      //                             Strain rate squared calculation
      //------------------------------------------------------------------------------------- Strain rate squared calculation
      strain_rate_tilde_square += strain_rate_tilde[i + DXQY::nD*j]*strain_rate_tilde[i + DXQY::nD*j];
    } 
  }

  rhoK_ = rhoKInput+lbBaseEps;
  rhoE_ = rhoEInput+lbBaseEps;
  

  const lbBase_t rhoInv = 1./rho;
  const lbBase_t gammaDotTildeSquare = 2*strain_rate_tilde_square;
  const lbBase_t Mk = DXQY::qSum(g);
  const lbBase_t Me = DXQY::qSum(h);


  //                             Solving for rhoK and rhoE
  //------------------------------------------------------------------------------------- Solving for rhoK and rhoE
  
  
  tau_ = rhoInv*X1_*rhoK_*rhoK_/rhoE_;

  sourceK_ = 0.0;
  sourceE_ = 0.0;
  

  const lbBase_t maxtau = 2.0; // 0.75  
  //if(tau_<maxtau){

  //                                       Alt. 1
  //------------------------------------------------------------------------------------- Alt 1.  
  const lbBase_t tauTauInvSquare = tau_/((tau0_ + tau_)*(tau0_ + tau_));
  sourceK_ = rhoInv*Y1_*gammaDotTildeSquare*tauTauInvSquare - rhoE_;
  sourceE_ = (rhoE_/rhoK_)*(Z1_*rhoInv*gammaDotTildeSquare*tauTauInvSquare - Z2_*rhoE_);
  
  

  /*
  //                                       Alt. 2
  //------------------------------------------------------------------------------------- Alt 1.  
  sourceK_ = 2*lbBaseEps;
  sourceE_ = 2*lbBaseEps;
  
  
  for (int i=0; i < 10; ++i){
    rhoK_ = Mk + 0.5*sourceK_;
    rhoE_ = Me + 0.5*sourceE_;
      
    tau_ = rhoInv*X1_*rhoK_*rhoK_/rhoE_;
   
    const lbBase_t tauTauInvSquare = tau_/((tau0_ + tau_)*(tau0_ + tau_));
    sourceK_ = rhoInv*Y1_*gammaDotTildeSquare*tauTauInvSquare - rhoE_;
    sourceE_ = (rhoE_/rhoK_)*(Z1_*rhoInv*gammaDotTildeSquare*tauTauInvSquare - Z2_*rhoE_);
  }
  */

  rhoK_ = Mk + 0.5*sourceK_;
  rhoE_ = Me + 0.5*sourceE_;

  if ( (rhoK_ < 0) || (rhoE_ < 0)) {
    sourceK_ = -2*Mk + 2*lbBaseEps;
    sourceE_ = -2*Me + 2*lbBaseEps;

    rhoK_ = Mk + 0.5*sourceK_;
    rhoE_ = Me + 0.5*sourceE_;
  }

  tau_ = rhoInv*X1_*rhoK_*rhoK_/rhoE_;
  if (tau_ < 0) {
    std::cout << tau_ << " " << rhoK_ << " " << rhoE_ << std::endl;
    MPI_Finalize();
    exit(1);
  }
  //}
  //else{
  //  tau_=maxtau;
  //}

  //                             Calculate relaxation times
  //------------------------------------------------------------------------------------- Calculate relaxation times 
  
  tauK_ = tauK0Term_ + 0.5 + tau_*sigmakInv_;
  tauE_ = tauE0Term_ + 0.5 + tau_*sigmaepsilonInv_;

  tau_ += tau0_;

  
  gammaDotTilde_= sqrt(2*strain_rate_tilde_square);
}



template <typename DXQY>
template <typename T1, typename T2, typename T3>
void Rans<DXQY>::zouHeFixedValueLeftBnd(
					const T1 &bndNodes,
					T2 &f,
					T3 &rho,
					const lbBase_t fixValue/*,
			                const T3 &grid*/ 		  
					)
/* Returns void: computes object values
 * 
 * 
 * 
 * Parameters
 * ----------
 * bndNodes : array-like [int], 
 *            list of boundary nodes
 *
 * f : array-like, size = [DXQY::nQ]
 *     lb distribution
 * 
 * fixValue : float-like
 *     scalar to be set at boundary
 * 
  
 *  : array-like, size = [DXQY::nD]
 *     body force
 * 
 * grid : grid object
 *    
 * 
 * Returns
 * -------
 * void
 *     
 */ 
{

  for (const auto &nodeNo: bndNodes) {
        const std::valarray<lbBase_t> fNode = f(0, nodeNo);
        //const std::valarray<lbBase_t> forceNode = force(0, nodeNo);
        const lbBase_t rho_ux = fixValue*rho(0, nodeNo) /*- 0.5*forceNode[0]*/ - (fNode[2] + fNode[6] + fNode[8] + 2*(fNode[3] + fNode[4] + fNode[5]));
        f(0, 0, nodeNo) = fNode[4] + (2./3.)*rho_ux /*- (1./3.)*forceNode[0]*/;
        f(0, 1, nodeNo) = fNode[5] + 0.5*(fNode[6] - fNode[2]) + (1./6.)*rho_ux /*+ (5./12.)*forceNode[0] + (1./4.)*forceNode[1]*/; 
        f(0, 7, nodeNo) = fNode[3] + 0.5*(fNode[2] - fNode[6]) + (1./6.)*rho_ux /*+ (5./12.)*forceNode[0] - (1./4.)*forceNode[1]*/;

    }
  
}

template <typename DXQY>
template <typename T1, typename T2/*, typename T3*/>
void Rans<DXQY>::zouHeFixedVelocityLeftBnd(
					   const T1 &bndNodes,
					   T2 &f,
					   const lbBase_t fixValue/*,
					   const T3 &grid*/ 		  
					   )
/* Returns void: computes object values
 * 
 * 
 * 
 * Parameters
 * ----------
 * bndNodes : array-like [int], 
 *            list of boundary nodes
 *
 * f : array-like, size = [DXQY::nQ]
 *     lb distribution
 * 
 * fixValue : float-like
 *     scalar to be set at boundary
 * 
  
 *  : array-like, size = [DXQY::nD]
 *     body force
 * 
 * grid : grid object
 *    
 * 
 * Returns
 * -------
 * void
 *     
 */ 
{

  for (const auto &nodeNo: bndNodes) {
        const std::valarray<lbBase_t> fNode = f(0, nodeNo);
	lbBase_t velNoise = velNoiseAmplitude_*((2.0*std::rand())/lbBase_t(RAND_MAX) - 1);
        //const std::valarray<lbBase_t> forceNode = force(0, nodeNo);	
        //const lbBase_t rho = 1/(1-(fixValue))*(fNode[2] + fNode[6] + fNode[8] + 2*(fNode[3] + fNode[4] + fNode[5]));
	const lbBase_t rho =1.0;
	const lbBase_t rho_ux = rho*fixValue+velNoise;
        f(0, 0, nodeNo) = fNode[4] + (2./3.)*rho_ux /*- (1./3.)*forceNode[0]*/;
        f(0, 1, nodeNo) = fNode[5] + 0.5*(fNode[6] - fNode[2]) + (1./6.)*rho_ux /*+ (5./12.)*forceNode[0] + (1./4.)*forceNode[1]*/; 
        f(0, 7, nodeNo) = fNode[3] + 0.5*(fNode[2] - fNode[6]) + (1./6.)*rho_ux /*+ (5./12.)*forceNode[0] - (1./4.)*forceNode[1]*/;

    }
  
}

template <typename DXQY>
template <typename T1, typename T2, typename T3, typename T4, typename T5>
void Rans<DXQY>::zouHeOpenRightBnd(
				   const T1 &bndNodes,
				   T2 &f,
				   //const lbBase_t fixValue,
				   const T3 &scalarField,
				   const T4 &force,
				   const T5 &grid 		  
				   )
/* Returns void: computes object values
 * 
 * 
 * 
 * Parameters
 * ----------
 * bndNodes : array-like [int], 
 *            list of boundary nodes
 *
 * f : array-like, size = [DXQY::nQ]
 *     lb distribution
 * 
 * fixValue : float-like
 *     scalar to be set at boundary
 * 
  
 *  : array-like, size = [DXQY::nD]
 *     body force
 * 
 * grid : grid object
 *    
 * 
 * Returns
 * -------
 * void
 *     
 */ 
{


  for (const auto &nodeNo : bndNodes)
    {
      const std::valarray<lbBase_t> fNode = f(0, nodeNo);
      const std::valarray<lbBase_t> forceNode = force(0, nodeNo);  
      lbBase_t scalar0Node = 2*scalarField(0, grid.neighbor(4, nodeNo))  - scalarField(0, grid.neighbor(4, grid.neighbor(4, nodeNo)));
      
      const lbBase_t rho_ux = scalar0Node - 0.5 * forceNode[0] - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[0] + fNode[1] + fNode[7]));
      f(0, 4, nodeNo) = fNode[0] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
      f(0, 3, nodeNo) = fNode[7] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] + (1. / 4.) * forceNode[1];
      f(0, 5, nodeNo) = fNode[1] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] - (1. / 4.) * forceNode[1];
      
      
    }

  
  
}

template <typename DXQY>
template <typename T1, typename T2, typename T3, typename T4>
void Rans<DXQY>::zouHeFixedValueRightBnd(
				   const T1 &bndNodes,
				   T2 &f,
				   const lbBase_t fixValue,
				   const T3 &force,
				   const T4 &grid 		  
				   )
/* Returns void: computes object values
 * 
 * 
 * 
 * Parameters
 * ----------
 * bndNodes : array-like [int], 
 *            list of boundary nodes
 *
 * f : array-like, size = [DXQY::nQ]
 *     lb distribution
 * 
 * fixValue : float-like
 *     scalar to be set at boundary
 * 
  
 *  : array-like, size = [DXQY::nD]
 *     body force
 * 
 * grid : grid object
 *    
 * 
 * Returns
 * -------
 * void
 *     
 */ 
{


  for (const auto &nodeNo : bndNodes)
    {
      const std::valarray<lbBase_t> fNode = f(0, nodeNo);
      const std::valarray<lbBase_t> forceNode = force(0, nodeNo);  
      lbBase_t scalar0Node = fixValue;
      
      const lbBase_t rho_ux = scalar0Node - 0.5 * forceNode[0] - (fNode[2] + fNode[6] + fNode[8] + 2 * (fNode[0] + fNode[1] + fNode[7]));
      f(0, 4, nodeNo) = fNode[0] + (2. / 3.) * rho_ux - (1. / 3.) * forceNode[0];
      f(0, 3, nodeNo) = fNode[7] + 0.5 * (fNode[6] - fNode[2]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] + (1. / 4.) * forceNode[1];
      f(0, 5, nodeNo) = fNode[1] + 0.5 * (fNode[2] - fNode[6]) + (1. / 6.) * rho_ux + (5. / 12.) * forceNode[0] - (1. / 4.) * forceNode[1];
      
      
    }

  
  
}

template <typename DXQY>
template <typename T1, typename T2>
void Rans<DXQY>::copyDistBnd(
				   const T1 &bndNodes,
				   T2 &f,
           int qDir,
           const Grid<DXQY> &grid 
           )
{
  for (const auto &nodeNo : bndNodes)
  {
    f.set(0, nodeNo) = f(0, grid.neighbor(qDir, nodeNo));
  }
}

template <typename DXQY>
lbBase_t Rans<DXQY>::wallShear(lbBase_t uIntrp)
{
  lbBase_t ret = 0;
  if (uIntrp < uWallMin_) {
    ret = std::sqrt(viscosity0_*uIntrp/yp_);
  } else if (uIntrp < uWallMax_) {
    ret = viscosity0_*155.0/yp_; // Binary search
    for (int n = 0; n < 4; ++n) {
      const lbBase_t AlnBx = lowC0_*std::log(lowC1_*ret);
      ret += -(ret*AlnBx - uIntrp)/(AlnBx - lowC0_);
    }
  } else {
    ret = viscosity0_*300.0/yp_;
    std::cout << "uInterp (> uMax) = " << uIntrp << std::endl;
  }

  return ret;
}

template <typename DXQY>
void Rans<DXQY>::solidBnd(
      LbField<DXQY> &f,
      LbField<DXQY> &fTmp,
      ScalarField &rho,
      VectorField<DXQY> &vel,
      ScalarField &viscocity,
      LbField<DXQY> &g,
      LbField<DXQY> &gTmp,
      ScalarField &rhoK, 
      LbField<DXQY> &h,
      LbField<DXQY> &hTmp,
      ScalarField &rhoE,
      std::vector<InterpolationElement> &boundarynodes, 
      Nodes<DXQY> &nodes,
      Grid<DXQY> &grid)
{
    int cnt = 0;
    lbBase_t sumrho = 0;
    for (auto &bn: boundarynodes) {
      sumrho += interpolateScalar(rho, bn.wb, bn.pnts);
      cnt += 1;
    }

    for (auto &bn: boundarynodes) {
        const int nodeNo = bn.nodeNo;
        const lbBase_t rhoNode = sumrho/cnt;//1.0; //interpolateScalar(rho, bn.wb, bn.pnts);

        std::valarray<lbBase_t> nvec = {bn.normal[0], bn.normal[1]};
        std::valarray<lbBase_t> tvec = {nvec[1], -nvec[0]};
        if (tvec[0] < 0) {
          tvec[0] = -tvec[0];
          tvec[1] = -tvec[1];
        }

        const std::valarray<lbBase_t> cn = DXQY::cDotAll(nvec);
        const std::valarray<lbBase_t> ct = DXQY::cDotAll(tvec);

        std::valarray<lbBase_t> shearStress_pnts(0.0, bn.pnts.size());
        lbBase_t shearStress_mean = 0;
        lbBase_t nuPnts_mean = 0;
        for (std::size_t pntNo = 0; pntNo < bn.pnts.size(); ++pntNo) {
          const int nodePntNo = bn.pnts[pntNo];
          const lbBase_t rhoPnt = rho(0, nodePntNo);
          const std::valarray<lbBase_t> velPnt = vel(0, nodePntNo);
          const std::valarray<lbBase_t> fNeqPnt = f(0, nodePntNo) - calcfeq<DXQY>(rhoPnt, velPnt);
          lbBase_t Mnt = 0.0;          
          for (int q=0; q < DXQY::nQ; ++q) {
            Mnt += cn[q]*ct[q]*fNeqPnt[q];
          }
          const lbBase_t nu = viscocity(0, nodePntNo);
          shearStress_pnts[pntNo] = -DXQY::c2Inv*nu*Mnt/(nu + 0.5);
          shearStress_mean +=  -bn.wb[pntNo]*shearStress_pnts[pntNo];
          nuPnts_mean += bn.wb[pntNo]*nu;
        }

        // Calculate y_plus
        lbBase_t rhoPnts_mean = interpolateScalar(rho, bn.wb, bn.pnts);
        lbBase_t dx_dut = 0.5*( DXQY::dot(vel(0, bn.pnts[1]), tvec) - DXQY::dot(vel(0, bn.pnts[0]), tvec) +
                                DXQY::dot(vel(0, bn.pnts[2]), tvec) - DXQY::dot(vel(0, bn.pnts[3]), tvec));
        lbBase_t dy_dut = 0.5*( DXQY::dot(vel(0, bn.pnts[3]), tvec) - DXQY::dot(vel(0, bn.pnts[0]), tvec) +
                                DXQY::dot(vel(0, bn.pnts[2]), tvec) - DXQY::dot(vel(0, bn.pnts[1]), tvec));

        shearStress_mean = rhoPnts_mean*nuPnts_mean*(nvec[0]*dx_dut + nvec[1]*dy_dut);

        lbBase_t u_star = std::sqrt(std::abs(shearStress_mean)/rhoPnts_mean);

        viscocity(0, nodeNo) = u_star; 
    }



    for (auto &bn: boundarynodes) {
        const int nodeNo = bn.nodeNo;
        const lbBase_t rhoNode = sumrho/cnt;//1.0; //interpolateScalar(rho, bn.wb, bn.pnts);

        std::valarray<lbBase_t> nvec = {bn.normal[0], bn.normal[1]};
        std::valarray<lbBase_t> tvec = {nvec[1], -nvec[0]};
        if (tvec[0] < 0) {
          tvec[0] = -tvec[0];
          tvec[1] = -tvec[1];
        }

        const std::valarray<lbBase_t> cn = DXQY::cDotAll(nvec);
        const std::valarray<lbBase_t> ct = DXQY::cDotAll(tvec);

        std::valarray<lbBase_t> shearStress_pnts(0.0, bn.pnts.size());
        lbBase_t shearStress_mean = 0;
        lbBase_t nuPnts_mean = 0;
        for (std::size_t pntNo = 0; pntNo < bn.pnts.size(); ++pntNo) {
          const int nodePntNo = bn.pnts[pntNo];
          const lbBase_t rhoPnt = rho(0, nodePntNo);
          const std::valarray<lbBase_t> velPnt = vel(0, nodePntNo);
          const std::valarray<lbBase_t> fNeqPnt = f(0, nodePntNo) - calcfeq<DXQY>(rhoPnt, velPnt);
          lbBase_t Mnt = 0.0;          
          for (int q=0; q < DXQY::nQ; ++q) {
            Mnt += cn[q]*ct[q]*fNeqPnt[q];
          }
          const lbBase_t nu = viscocity(0, nodePntNo);
          shearStress_pnts[pntNo] = -DXQY::c2Inv*nu*Mnt/(nu + 0.5);
          shearStress_mean +=  -bn.wb[pntNo]*shearStress_pnts[pntNo];
          nuPnts_mean += bn.wb[pntNo]*nu;
        }

        // Calculate y_plus
        lbBase_t rhoPnts_mean = interpolateScalar(rho, bn.wb, bn.pnts);
        lbBase_t dx_dut = 0.5*( DXQY::dot(vel(0, bn.pnts[1]), tvec) - DXQY::dot(vel(0, bn.pnts[0]), tvec) +
                                DXQY::dot(vel(0, bn.pnts[2]), tvec) - DXQY::dot(vel(0, bn.pnts[3]), tvec));
        lbBase_t dy_dut = 0.5*( DXQY::dot(vel(0, bn.pnts[3]), tvec) - DXQY::dot(vel(0, bn.pnts[0]), tvec) +
                                DXQY::dot(vel(0, bn.pnts[2]), tvec) - DXQY::dot(vel(0, bn.pnts[1]), tvec));

        shearStress_mean = rhoPnts_mean*nuPnts_mean*(nvec[0]*dx_dut + nvec[1]*dy_dut);

        lbBase_t u_star = std::sqrt(std::abs(shearStress_mean)/rhoPnts_mean);

        lbBase_t numNeig = 0;
        u_star = 0;
        for (auto neigNo:grid.neighbor(nodeNo)) {
          if (nodes.getType(neigNo) == 2) {
            numNeig += 1.0;
            u_star += viscocity(0, neigNo);
          }
        }
        u_star /= numNeig;

        lbBase_t y_plus = yp_*u_star/viscosity0_;
        // Calculate velocity at the wall
        lbBase_t u_wall = 0;
        if (y_plus < 10.92) {
          u_wall = u_star*y_plus;
        } 
        else if (y_plus < 300 ) {
          u_wall = u_star*std::log(E_*y_plus)/kappa_;
        }
        else 
        {
          std::cout << "Warning y_plus = " << y_plus << std::endl;
          u_wall = std::log(E_*300)/kappa_;
        }
        // Calculate velocity at the boundary
        const std::valarray<lbBase_t> velPnts_mean = interpolateVector(vel, bn.wb, bn.pnts); // Bulk fluid
        std::valarray<lbBase_t> uNode(0.0, DXQY::nD);
        if (u_wall > std::abs(DXQY::dot(tvec, velPnts_mean)))  u_wall = std::abs(DXQY::dot(tvec, velPnts_mean));
        if (DXQY::dot(tvec, velPnts_mean) > 0) {
          uNode = (bn.gamma*velPnts_mean + bn.gamma2*u_wall*tvec)/(bn.gamma + bn.gamma2);
        } 
        else {
          uNode = (bn.gamma*velPnts_mean - bn.gamma2*u_wall*tvec)/(bn.gamma + bn.gamma2);
        } 

        const lbBase_t k_wall = u_star*u_star/std::sqrt(Cmu_);
        const lbBase_t epsilon_wall = u_star*u_star*u_star/(kappa_*yp_);

        // turbulent kinematic viscosity at yp
        //const lbBase_t nu_t_wall = (epsilon_wall > 0) ? Cmu_*k_wall*k_wall/epsilon_wall : 0.0;
        const lbBase_t nu_t_wall = kappa_*yp_*u_star;
        const lbBase_t nu_wall = nu_t_wall + viscosity0_;
        const lbBase_t nuNode = (bn.gamma*nuPnts_mean + bn.gamma2*nu_wall)/(bn.gamma + bn.gamma2);

        const lbBase_t const_fneq = DXQY::c2Inv*shearStress_mean*(DXQY::c2Inv*nuNode + 0.5)/nuNode;

        std::valarray<lbBase_t> fneq_wall(0.0, DXQY::nQ);
        for (int q=0; q < DXQY::nQ; ++q) {
          fneq_wall[q] = -DXQY::w[q]*cn[q]*ct[q]*const_fneq;
        }

        const std::valarray<lbBase_t> fneqNode = (bn.gamma*interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts) + bn.gamma2*fneq_wall)/(bn.gamma + bn.gamma2); 
        // interpolated velocity
        // equilibrium distribution
        const std::valarray<lbBase_t> feqNode =  calcfeq<DXQY>(rhoNode, uNode);

        /*f.set(0, nodeNo) = feqNode + 0*fneqNode;

        const lbBase_t rhoKNode =  (bn.gamma*interpolateScalar(rhoK, bn.wb, bn.pnts) + bn.gamma2*rhoNode*k_wall)/(bn.gamma + bn.gamma2); 

        g.set(0, nodeNo) = calcfeq<DXQY>(rhoKNode, uNode) + interpolateLBFieldNeq(rhoK, vel, g, bn.wb, bn.pnts);

        const lbBase_t rhoEpsilonNode = (bn.gamma*interpolateScalar(rhoE, bn.wb, bn.pnts) + bn.gamma2*rhoNode*epsilon_wall)/(bn.gamma + bn.gamma2);

        h.set(0, nodeNo) = calcfeq<DXQY>(rhoEpsilonNode, uNode) + interpolateLBFieldNeq(rhoE, vel, h, bn.wb, bn.pnts);
        */

        const std::valarray<lbBase_t> fNode = feqNode + 0*fneqNode;
        const lbBase_t tauNode = DXQY::c2Inv*nuNode + 0.5;
        // Run one iteration 
        const lbBase_t u2 = DXQY::dot(uNode, uNode);
        const auto cu = DXQY::cDotAll(uNode);
        const auto omegaBGK = calcOmegaBGKTRT<DXQY>(fNode, tauNode, 1.0, rhoNode, u2, cu);
        fTmp.propagateTo(0, nodeNo, fNode + omegaBGK, grid);
        rho(0, nodeNo) = rhoNode;
        vel.set(0, nodeNo) = uNode;


        const lbBase_t rhoKNode =  (bn.gamma*interpolateScalar(rhoK, bn.wb, bn.pnts) + bn.gamma2*rhoNode*k_wall)/(bn.gamma + bn.gamma2);
        const lbBase_t tauKNode = tauK0Term_ + 0.5 + (tauNode - 0.5 - tau0_)*sigmakInv_;
        const std::valarray<lbBase_t> gNode = calcfeq<DXQY>(rhoKNode, uNode) + 0* interpolateLBFieldNeq(rhoK, vel, g, bn.wb, bn.pnts);
        const auto omegaBGK_K = calcOmegaBGKTRT<DXQY>(gNode, 1.0, tauKNode, rhoKNode, u2, cu);

        const lbBase_t rhoENode = (bn.gamma*interpolateScalar(rhoE, bn.wb, bn.pnts) + bn.gamma2*rhoNode*epsilon_wall)/(bn.gamma + bn.gamma2);
        const lbBase_t tauENode = tauE0Term_ + 0.5 + (tauNode - 0.5 - tau0_)*sigmaepsilonInv_;
        const std::valarray<lbBase_t> hNode  = calcfeq<DXQY>(rhoENode, uNode) + 0* interpolateLBFieldNeq(rhoE, vel, h, bn.wb, bn.pnts); 
        const auto omegaBGK_E = calcOmegaBGKTRT<DXQY>(hNode, 1.0, tauENode, rhoENode, u2, cu);

      //                               Collision and propagation
      //------------------------------------------------------------------------------------- Collision and propagation
      fTmp.propagateTo(0, nodeNo, fNode + omegaBGK, grid);
      gTmp.propagateTo(0, nodeNo, gNode + omegaBGK_K, grid);
      hTmp.propagateTo(0, nodeNo, hNode + omegaBGK_E, grid);



    }
}


















template <typename DXQY>
void Rans<DXQY>::solidBndVel(
      LbField<DXQY> &f,
      ScalarField &rho,
      VectorField<DXQY> &vel,
      ScalarField &viscocity,
      LbField<DXQY> &g,
      ScalarField &rhoK, 
      LbField<DXQY> &h,
      ScalarField &rhoE,
      std::vector<InterpolationElement> &boundarynodes, 
      Grid<DXQY> &grid)
{

    for (auto &bn: boundarynodes) {
        const int nodeNo = bn.nodeNo;
        const lbBase_t rhoNode = 1.0; //interpolateScalar(rho, bn.wb, bn.pnts);
        // velNodeInterp: extrapolated value at yp
        // tvec: tangent to the surface
        const std::valarray<lbBase_t> tvec = {-bn.normal[1], bn.normal[0]};
        const std::valarray<lbBase_t> velNodeIntrpb = interpolateVector(vel, bn.wb, bn.pnts); // Bulk fluid
        const std::valarray<lbBase_t> velNodeIntrpa = interpolateVector(vel, bn.wa, bn.pnts); // Wall
        const std::valarray<lbBase_t> velNodeIntrp = std::abs(DXQY::dot(velNodeIntrpa, tvec)) < std::abs(DXQY::dot(velNodeIntrpb, tvec)) ? velNodeIntrpa : velNodeIntrpb; 

        // Tangent of the extrapolated value
        lbBase_t velWall_t = DXQY::dot(velNodeIntrp, tvec);

        // If the interpolated value changes sign we will just set it to zero
        if (velWall_t*DXQY::dot(velNodeIntrpb, tvec) < 0) {
                velWall_t = 0;
        }  
        // Assumes that the velocity should decrease (or at least not increase) as you go towards 
        // the wall. (Need to check if this the case in transient flows)
        if (  std::abs(velWall_t) >  std::abs(DXQY::dot(velNodeIntrp, tvec)) ) {
            velWall_t = DXQY::dot(velNodeIntrp, tvec);
        }

        velWall_t = 0;

        // velWall is the value that is used for yp (that is at the wall in our case)
        std::valarray<lbBase_t> velWall_y = velWall_t*tvec;  
        const lbBase_t uStar = wallShear(std::abs(velWall_t));
        // wallShearStress: (Shear stress at yp but also at the actual wall)
        const lbBase_t wallShearStress = rhoNode*uStar*uStar;
        // Diffusive fields at yp
        const lbBase_t k_y = uStar*uStar/std::sqrt(Cmu_);
        const lbBase_t epsilon_y = uStar*uStar*uStar/(kappa_*yp_);

        // turbulent kinematic viscosity at yp
        const lbBase_t my_t_y = (epsilon_y > 0) ? Cmu_*k_y*k_y/epsilon_y : 0.0;
        // Calculate f_neq
        // -- relaxation time - 0.5 at yp
        lbBase_t  tmp_fneq = DXQY::c2Inv*(my_t_y + viscosity0_);
        // -- calculate Q_nt = (c*n)(c*t)
        const std::valarray<lbBase_t> cn = DXQY::cDotAll(bn.normal);
        const std::valarray<lbBase_t> ct = DXQY::cDotAll(tvec);
        std::valarray<lbBase_t> fNeq_y(DXQY::nQ);
        // -- tmp_fneq = rho*tau*tau_ns/(tau-0.5)/cs_4
        tmp_fneq = DXQY::c4Inv*rhoNode*(tmp_fneq + 0.5)*wallShearStress/tmp_fneq;
        if (velWall_t > 0) tmp_fneq = -tmp_fneq;
        for (int q=0; q < DXQY::nQ; ++q) {
          fNeq_y[q] = -DXQY::w[q]*cn[q]*ct[q]*tmp_fneq; 
        }
/*        const std::valarray<lbBase_t> fneqTmp = interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts);
        lbBase_t tmp00 = 0;
        lbBase_t tmp10 = 0;
        lbBase_t tmp01 = 0;
        lbBase_t tmp11 = 0;

        lbBase_t amp00 = 0;
        lbBase_t amp10 = 0;
        lbBase_t amp01 = 0;
        lbBase_t amp11 = 0;


        for (int q=0; q < DXQY::nQ; ++q) {
          lbBase_t c0 =  D2Q9::c(q, 0);
          lbBase_t c1 =  D2Q9::c(q, 1);
          tmp00 += (c0*c0 - 1.0/3.0)*fNeq_y[q];
          tmp10 += c1*c0*fNeq_y[q];
          tmp11 += (c1*c1 - 1.0/3.0)*fNeq_y[q];
          amp00 += (c0*c0 - 1.0/3.0)*fneqTmp[q];
          amp10 += c1*c0*fneqTmp[q];
          amp11 += (c1*c1 - 1.0/3.0)*fneqTmp[q];
        }
        lbBase_t n0 = bn.normal[0];
        lbBase_t n1 = bn.normal[1];
        lbBase_t t0 = tvec[0];
        lbBase_t t1 = tvec[1];
        lbBase_t tcalc = tmp00*n0*t0 + tmp01*(n0*t1 + n1*t0) + tmp11*n1*t1;
        lbBase_t acalc = amp00*n0*t0 + amp01*(n0*t1 + n1*t0) + amp11*n1*t1;
        std::cout << tcalc << " " << acalc << " " <<  wallShearStress << " " << t0 << std::endl; */

        // interpolated non-equilibrium distribution
        const std::valarray<lbBase_t> fneqNode = (bn.gamma*interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts) + bn.gamma2*fNeq_y)/(bn.gamma + bn.gamma2); 
        // interpolated velocity
        const std::valarray<lbBase_t> velNode = (bn.gamma*interpolateVector(vel, bn.wb, bn.pnts) + bn.gamma2*velWall_y)/(bn.gamma + bn.gamma2);
        // equilibrium distribution
        const std::valarray<lbBase_t> feqNode =  calcfeq<DXQY>(rhoNode, velNode);

        f.set(0, nodeNo) = feqNode + fneqNode;

        const lbBase_t rhoKNode =  (bn.gamma*interpolateScalar(rhoK, bn.wb, bn.pnts) + bn.gamma2*rhoNode*k_y)/(bn.gamma + bn.gamma2); 

        g.set(0, nodeNo) = calcfeq<DXQY>(rhoKNode, velNode) + interpolateLBFieldNeq(rhoK, vel, g, bn.wb, bn.pnts);;

        const lbBase_t rhoEpsilonNode = (bn.gamma*interpolateScalar(rhoE, bn.wb, bn.pnts) + bn.gamma2*rhoNode*epsilon_y)/(bn.gamma + bn.gamma2);;

        h.set(0, nodeNo) = calcfeq<DXQY>(rhoEpsilonNode, velNode) + interpolateLBFieldNeq(rhoE, vel, h, bn.wb, bn.pnts); ;

        // std::cout << k_y << "(" <<  rhoKNode <<")" << " " << epsilon_y << "(" <<  rhoEpsilonNode <<")"  << "  " << yp_*rhoNode*uStar/viscosity0_  << std::endl;

        /* const std::valarray<lbBase_t> velNode = (bn.gamma*velNodeIntrp + bn.gamma2*velWall)/(bn.gamma + bn.gamma2);

        //const std::valarray<lbBase_t> velNode = bn.gamma*velNodeIntrp/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> fneq = interpolateLBFieldNeq(rho, vel, f, bn.wb, bn.pnts);
        const lbBase_t dvel_dn = DXQY::dot(velNodeIntrp, tvec)/(bn.gamma + bn.gamma2);
        const std::valarray<lbBase_t> forceNode = rampup*bodyForce(0, nodeNo);
        const std::valarray<lbBase_t> fneq_wall = calcfneqWall<DXQY>(tau, rhoNode, dvel_dn, forceNode, bn);

        const std::valarray<lbBase_t> delta_fneq = bn.gamma2*calcDeltafneqwall<DXQY>(fneq, bn)/(bn.gamma + bn.gamma2);

        f.set(0, nodeNo) = calcfeq<DXQY>(rhoNode, velNode) - 3*0.5*D2Q9::cDotAll(forceNode) + 0*(bn.gamma*fneq + bn.gamma2*fneq_wall)/(bn.gamma + bn.gamma2);
        */
    }



}


#endif

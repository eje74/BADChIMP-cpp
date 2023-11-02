#ifndef LBRANS_H
#define LBRANS_H

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../lbsolver/LBglobal.h"
#include "../io/Input.h"

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
  

private:
  //                             Constant input parameters
  //------------------------------------------------------------------------------------- Constant input parameters
    const lbBase_t tau0_;
    const lbBase_t Cmu_;
    const lbBase_t C1epsilon_;
    const lbBase_t C2epsilon_;
    const lbBase_t sigma0kInv_;
    const lbBase_t sigmakInv_;
    const lbBase_t sigma0epsilonInv_;
    const lbBase_t sigmaepsilonInv_;

    const lbBase_t velNoiseAmplitude_;
  //------------------------------------------------------------------------------------- Constant derived varaibles
    const lbBase_t X1_;
    const lbBase_t Y1_;
    const lbBase_t Z1_;
    const lbBase_t Z2_;

    const lbBase_t tauK0Term_;
    const lbBase_t tauE0Term_;
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
  :tau0_(DXQY::c2Inv * input["fluid"]["viscosity"] + 0.5), 
   Cmu_(input["RANS"]["k-epsilonCoef"]["C_mu"]),
   C1epsilon_(input["RANS"]["k-epsilonCoef"]["C_1epsilon"]),
   C2epsilon_(input["RANS"]["k-epsilonCoef"]["C_2epsilon"]),
   sigma0kInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_0k"]),
   sigmakInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_k"]),
   sigma0epsilonInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_0epsilon"]),
   sigmaepsilonInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_epsilon"]),
   velNoiseAmplitude_(input["RANS"]["inlet"]["velNoiseAmplitude"]),
   X1_(Cmu_*DXQY::c2Inv),
   Y1_(0.25*DXQY::c2Inv),
   Z1_(0.25*DXQY::c2Inv*C1epsilon_),
   Z2_(C2epsilon_),
   tauK0Term_((tau0_-0.5)*sigma0kInv_),
   tauE0Term_((tau0_-0.5)*sigma0epsilonInv_)

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
  std::cout << "X1 = " << X1_ << std::endl;
  std::cout << "Y1 = " << Y1_ << std::endl;
  std::cout << "Z1 = " << Z1_ << std::endl;
  std::cout << "Z2 = " << Z2_ << std::endl;
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
  
  if(tau_<0.75){

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

  tau_ = rhoInv*X1_*rhoK_*rhoK_/rhoE_;
  }
  else{
    tau_=0.75;
  }

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



#endif

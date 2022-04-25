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
    Rans(Input &input);

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
	     const T1 &h);

    inline lbBase_t tau() const {return tau_;}
    inline lbBase_t tauK() const {return tauK_;}
    inline lbBase_t tauE() const {return tauE_;}

    inline lbBase_t sourceK() const {return sourceK_;}
    inline lbBase_t sourceE() const {return sourceE_;}

    inline lbBase_t rhoK() const {return rhoK_;}
    inline lbBase_t rhoE() const {return rhoE_;}

  inline lbBase_t gammaDot() const {return gammaDotTilde_;}

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
Rans<DXQY>::Rans(Input &input)
  :tau0_(DXQY::c2Inv * input["fluid"]["viscosity"] + 0.5), 
   Cmu_(input["RANS"]["k-epsilonCoef"]["C_mu"]),
   C1epsilon_(input["RANS"]["k-epsilonCoef"]["C_1epsilon"]),
   C2epsilon_(input["RANS"]["k-epsilonCoef"]["C_2epsilon"]),
   sigma0kInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_0k"]),
   sigmakInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_k"]),
   sigma0epsilonInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_0epsilon"]),
   sigmaepsilonInv_(1./input["RANS"]["k-epsilonCoef"]["sigma_epsilon"]),
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
		       const T1 &h)
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

  
  //gammaDotTilde_= sqrt(2*strain_rate_tilde_square);

  const lbBase_t rhoInv = 1./rho;
  const lbBase_t gammaDotTildeSquare = 2*strain_rate_tilde_square;
  lbBase_t tauTauInvSquare;
  rhoK_ = DXQY::qSum(g)+ lbBaseEps;
  rhoE_ = DXQY::qSum(h) + lbBaseEps;

  tau_ = rhoInv*X1_*rhoK_*rhoK_/rhoE_;

  for (int i=0; i < 10; ++i){
    tauTauInvSquare = tau_/((tau0_ + tau_)*(tau0_ + tau_));
    sourceK_ = rhoInv*Y1_*gammaDotTildeSquare*tauTauInvSquare - rhoE_;
    rhoK_ += 0.5*sourceK_;
    sourceE_ = rhoInv*(rhoE_/rhoK_)*(Z1_*gammaDotTildeSquare*tauTauInvSquare - Z2_*rhoE_);
    rhoE_ += 0.5*sourceE_;
    tau_ = rhoInv*X1_*rhoK_*rhoK_/rhoE_;  // tau_ = \tau_t 
  }
  
  //                             Calculate relaxation times
  //------------------------------------------------------------------------------------- Calculate relaxation times 

  tauK_ = tauK0Term_ + 0.5 + tau_*sigmakInv_;
  tauE_ = tauE0Term_ + 0.5 + tau_*sigmaepsilonInv_;

  tau_ += tau0_;
  

  
}

#endif

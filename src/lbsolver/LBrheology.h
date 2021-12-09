#ifndef LBPOLYMER_H
#define LBPOLYMER_H

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../lbsolver/LBglobal.h"


template <typename DXQY>
class Newtonian
/* Class for Newtonian fluids
 * 
 * Attributes
 * ----------
 * 
 * Methods
 * -------
 * tau()
 *     returns the relaxation time calculated in ``omegaBGK``
 *
 * viscosity()
 *     returns the viscosity calculated in ``omegaBGK``
 *
 * EE()
 *     returns the strainrate constration calculated in "omegaBGK"
 *
 * omegaBGK(f, rho, u, u_sq, cu, force, source)
 *     returns the BKG collision term 
 */ 
{
public:
    Newtonian(lbBase_t tauInit);
    template <typename T1, typename T2, typename T3>
    std::valarray<lbBase_t> omegaBGK(const T1 &f,
                                     const lbBase_t& rho,
                                     const T2 &u,
                                     const lbBase_t& u_sq,
                                     const std::valarray<lbBase_t> &cu,
                                     const T3 &force,
                                     const lbBase_t &source);
    inline lbBase_t tau() const {
        return tau_;
    }
    inline lbBase_t viscosity() const {
        return visc_;
    }
    inline lbBase_t gammaDot() const {
        return gammaDot_;
    }
    
private:
    lbBase_t tau_;
    lbBase_t visc_;
    lbBase_t gammaDot_;
    std::vector<lbBase_t> tabular_strain_rate_;
};


template <typename DXQY>
Newtonian<DXQY>::Newtonian(lbBase_t tauInit)
/* Class constructor, sets object tau value from input tau
 * 
 * Here assume that we are given the viscosity as a function of a strain rate like parameter:
 *    \gamma = \sqrt{E_{ij}E_{ij}},
 * where
 *    E_{ij} = -\left[\sum_\alpha f_\alpha^\mathrm{neq} Q_{\alpha ij} + \frac{1}{2} \big(u_iF_j+u_jF_i + u_i u_j q \big)\Delta t\right]
 * so that
 *    tau = = 1/(\rho c^2_\mathrm{s} \Delta t)\mu_eff(\gamma) + 1/2
 * 
 * Parameters
 * ----------
 * tauInit : lbBase_t initial tau value given in input file.
 * 
 * Returns
 * -------
 * Newtonian<DXQY> constructor
 *    DXQY is a lattice type
 */ 
{
  tau_ = tauInit;
}


template <typename DXQY>
template <typename T1, typename T2, typename T3>
std::valarray<lbBase_t> Newtonian<DXQY>::omegaBGK(
                                 const T1 &f,
                                 const lbBase_t& rho,
                                 const T2 &u,
                                 const lbBase_t& u_sq,
                                 const std::valarray<lbBase_t> &cu,
                                 const T3 &F,
                                 const lbBase_t &source)
/* Returns the omegaBGK part of the lattice Boltzmann equation
 * 
 * This function calculates and returns the  collision operator is given by
 *     OmegaBGK = 1/tau*(f_\alpha - f_\alpha^neq),
 * where tau is the strain rate dependent relaxation time, f_\alpha^neq is the equilibirum  
 * consentration, f_\alpha is the lb distribution and \alpha is the basis velocity direction.
 *     The relaxation time tau and the viscosity are both calcualted in this functiona and
 * can be retrived using the class methods ``tau()`` and ``viscosity()``.
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
 * valarray, size = [DXQY::nQ]
 *     The BGK collision operator
 */ 
{
    std::valarray<lbBase_t> feq(DXQY::nQ);
    std::valarray<lbBase_t> strain_rate_tilde(0.0, DXQY::nD * DXQY::nD);
 
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        feq[q] = DXQY::w[q]*rho*(1 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq));        
        for (int i=0; i < DXQY::nD; ++i)
        {
            strain_rate_tilde[i + DXQY::nD*i] +=  (f[q] - feq[q])*(DXQY::c(q, i)*DXQY::c(q, i) - DXQY::c2);
            for (int j = i+1; j < DXQY::nD; ++j)
            {
                strain_rate_tilde[i + DXQY::nD*j] +=  (f[q] - feq[q])*DXQY::c(q, i)*DXQY::c(q, j);
            }            
        }
    } 
    
    lbBase_t strain_rate_tilde_square = 0;    
    for (int i=0; i < DXQY::nD; ++i)
    {
        strain_rate_tilde[i + DXQY::nD*i] += (u[i]*F[i]);
        strain_rate_tilde[i + DXQY::nD*i] += 0.5*u_sq * source;
        // Strain rate calculation
        strain_rate_tilde_square += strain_rate_tilde[i + DXQY::nD*i]*strain_rate_tilde[i + DXQY::nD*i];
        for (int j = i+1; j < DXQY::nD; ++j)
        {
            strain_rate_tilde[i + DXQY::nD*j] += 0.5*(u[i]*F[j] + u[j]*F[i]);
            strain_rate_tilde[i + DXQY::nD*j] += 0.5*u[i]*u[j]*source;
            // Strain rate calculation
            strain_rate_tilde_square += 2*strain_rate_tilde[i + DXQY::nD*j]*strain_rate_tilde[i + DXQY::nD*j];
        }
    } 

    lbBase_t strain_rate_tilde_cubed = 0;  

    
    

    auto tau_inv = 1.0/tau_;
    gammaDot_= sqrt(2*strain_rate_tilde_square)/(2*rho)*tau_inv*DXQY::c2Inv;
      
    std::valarray<lbBase_t> omega(DXQY::nQ);
    
    for (int q=0; q < DXQY::nQ; ++q)
    {
        omega[q] = -tau_inv*(f[q] - feq[q]); 
    }
        
    return omega;
}
//-------------------------------------------------------

template <typename DXQY>
class GeneralizedNewtonian
/* Class for generalized Newtonian fluids
 * 
 * Attributes
 * ----------
 * 
 * Methods
 * -------
 * tau()
 *     returns the relaxation time calculated in ``omegaBGK``
 *
 * viscosity()
 *     returns the viscosity calculated in ``omegaBGK``
 * 
 * omegaBGK(f, rho, u, u_sq, cu, force, source)
 *     returns the BKG collision term 
 */ 
{
public:
    GeneralizedNewtonian(std::string file_name);
    template <typename T1, typename T2, typename T3>
    std::valarray<lbBase_t> omegaBGK(const T1 &f,
                                     const lbBase_t& rho,
                                     const T2 &u,
                                     const lbBase_t& u_sq,
                                     const std::valarray<lbBase_t> &cu,
                                     const T3 &force,
                                     const lbBase_t &source);
    inline lbBase_t tau() const {
        return tau_;
    }
    inline lbBase_t viscosity() const {
        return visc_;
    }
    inline lbBase_t gammaDot() const {
        return gammaDot_;
    }
    void print_table() {
        std::cout << "SIZE = " << tabular_viscosity_.size() << std::endl;
        for (int i = 0; i < tabular_viscosity_.size(); ++i) {
            std::cout << tabular_strain_rate_[i] << " " << tabular_viscosity_[i] << std::endl;
        }
    }
private:
    lbBase_t tau_;
    lbBase_t visc_;
    lbBase_t gammaDot_;
    std::vector<lbBase_t> tabular_strain_rate_;
    std::vector<lbBase_t> tabular_viscosity_;
};


template <typename DXQY>
GeneralizedNewtonian<DXQY>::GeneralizedNewtonian(std::string file_name)
/* Class constructor, fills the share rate vs viscosity lookup-table
 * 
 * Here assume that we are given the viscosity as a function of a strain rate like parameter:
 *    \gamma = \sqrt{E_{ij}E_{ij}},
 * where
 *    E_{ij} = -\left[\sum_\alpha f_\alpha^\mathrm{neq} Q_{\alpha ij} + \frac{1}{2} \big(u_iF_j+u_jF_i + u_i u_j q \big)\Delta t\right]
 * so that
 *    tau = = 1/(\rho c^2_\mathrm{s} \Delta t)\mu_eff(\gamma) + 1/2
 * 
 * Parameters
 * ----------
 * file_name : string-like object
 *     filename including file path.
 * 
 * Returns
 * -------
 * GeneralizedNewtonian<DXQY> constructor
 *    DXQY is a lattice type
 */ 
{
    std::ifstream data_table_file(file_name, std::ios::in);

    int table_length;
    data_table_file >> table_length;

    tabular_strain_rate_.resize(table_length + 2);
    tabular_viscosity_.resize(table_length + 2);

    // Fill table
    for (int i = 0; i < table_length; ++i) {
        data_table_file >> tabular_strain_rate_[i+1] >> tabular_viscosity_[i+1];
    }
    // Lower bound
    tabular_strain_rate_[0] = tabular_strain_rate_[1] - 1;
    tabular_viscosity_[0] = tabular_viscosity_[1];
    // Upper bound
    tabular_strain_rate_[table_length + 1] = tabular_strain_rate_[table_length] + 1;
    tabular_viscosity_[table_length + 1] = tabular_viscosity_[table_length];
    // Close file
    data_table_file.close();
}


template <typename DXQY>
template <typename T1, typename T2, typename T3>
std::valarray<lbBase_t> GeneralizedNewtonian<DXQY>::omegaBGK(
                                 const T1 &f,
                                 const lbBase_t& rho,
                                 const T2 &u,
                                 const lbBase_t& u_sq,
                                 const std::valarray<lbBase_t> &cu,
                                 const T3 &F,
                                 const lbBase_t &source)
/* Returns the omegaBGK part of the lattice Boltzmann equation
 * 
 * This function calculates and returns the  collision operator is given by
 *     OmegaBGK = 1/tau*(f_\alpha - f_\alpha^neq),
 * where tau is the strain rate dependent relaxation time, f_\alpha^neq is the equilibirum  
 * consentration, f_\alpha is the lb distribution and \alpha is the basis velocity direction.
 *     The relaxation time tau and the viscosity are both calcualted in this functiona and
 * can be retrived using the class methods ``tau()`` and ``viscosity()``.
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
 * valarray, size = [DXQY::nQ]
 *     The BGK collision operator
 */ 
{
    std::valarray<lbBase_t> feq(DXQY::nQ);
    std::valarray<lbBase_t> strain_rate_tilde(0.0, DXQY::nD * DXQY::nD);
 
    for (int q = 0; q < DXQY::nQ; ++q)
    {
        feq[q] = DXQY::w[q]*rho*(1 + DXQY::c2Inv*cu[q] + DXQY::c4Inv0_5*(cu[q]*cu[q] - DXQY::c2*u_sq));        
        for (int i=0; i < DXQY::nD; ++i)
        {
            strain_rate_tilde[i + DXQY::nD*i] +=  (f[q] - feq[q])*(DXQY::c(q, i)*DXQY::c(q, i) - DXQY::c2);
            for (int j = i+1; j < DXQY::nD; ++j)
            {
                strain_rate_tilde[i + DXQY::nD*j] +=  (f[q] - feq[q])*DXQY::c(q, i)*DXQY::c(q, j);
            }            
        }
    } 
    
    lbBase_t strain_rate_tilde_square = 0;    
    for (int i=0; i < DXQY::nD; ++i)
    {
        strain_rate_tilde[i + DXQY::nD*i] += (u[i]*F[i]);
        strain_rate_tilde[i + DXQY::nD*i] += 0.5*u_sq * source;
        // Strain rate calculation
        strain_rate_tilde_square += strain_rate_tilde[i + DXQY::nD*i]*strain_rate_tilde[i + DXQY::nD*i];
        for (int j = i+1; j < DXQY::nD; ++j)
        {
            strain_rate_tilde[i + DXQY::nD*j] += 0.5*(u[i]*F[j] + u[j]*F[i]);
            strain_rate_tilde[i + DXQY::nD*j] += 0.5*u[i]*u[j]*source;
            // Strain rate calculation
            strain_rate_tilde_square += 2*strain_rate_tilde[i + DXQY::nD*j]*strain_rate_tilde[i + DXQY::nD*j];
        }
    } 
    
    // Lookup viscosity value
    auto upper = std::upper_bound(tabular_strain_rate_.begin()+1, tabular_strain_rate_.end()-1, strain_rate_tilde_square);
    auto i = std::distance(tabular_strain_rate_.begin()+1, upper);
    visc_ = (tabular_viscosity_[i+1] - tabular_viscosity_[i])*(strain_rate_tilde_square - tabular_strain_rate_[i])/(tabular_strain_rate_[i+1] - tabular_strain_rate_[i]) + tabular_viscosity_[i];
    tau_ = visc_*DXQY::c2Inv/rho + 0.5;

    auto tau_inv = 1.0/tau_;
    gammaDot_= sqrt(2*strain_rate_tilde_square)/(2*rho)*tau_inv*DXQY::c2Inv;
    
    std::valarray<lbBase_t> omega(DXQY::nQ);
    for (int q=0; q < DXQY::nQ; ++q)
    {
        omega[q] = -tau_inv*(f[q] - feq[q]); 
    }
        
    return omega;
}


#endif

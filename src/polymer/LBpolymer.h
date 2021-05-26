#ifndef LBPOLYMER_H
#define LBPOLYMER_H

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../lbsolver/LBglobal.h"


template <typename DXQY>
class GeneralizedNewtonian
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
    inline const lbBase_t tau() {
        return tau_;
    }
    inline const lbBase_t viscosity() {
        return visc_;
    }
    void print_table() {
        std::cout << "SIZE = " << viscosity_.size() << std::endl;
        for (int i = 0; i < viscosity_.size(); ++i) {
            std::cout << strain_rate_[i] << " " << viscosity_[i] << std::endl;
        }
    }
private:
    lbBase_t tau_;
    lbBase_t visc_;
    std::vector<lbBase_t> strain_rate_;
    std::vector<lbBase_t> viscosity_;
};


template <typename DXQY>
GeneralizedNewtonian<DXQY>::GeneralizedNewtonian(std::string file_name)
{
    std::ifstream data_table_file(file_name, std::ios::in);

    int table_length;
    data_table_file >> table_length;

    strain_rate_.resize(table_length + 2);
    viscosity_.resize(table_length + 2);

    // Fill table
    for (int i = 0; i < table_length; ++i) {
        data_table_file >> strain_rate_[i+1] >> viscosity_[i+1];
    }
    // Lower bound
    strain_rate_[0] = strain_rate_[1] - 1;
    viscosity_[0] = viscosity_[1];
    // Upper bound
    strain_rate_[table_length + 1] = strain_rate_[table_length] + 1;
    viscosity_[table_length + 1] = viscosity_[table_length];
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
 * OmegaBGK = 1/tau*(f_\alpha - f_\alpha^neq)
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
    auto upper = std::upper_bound(strain_rate_.begin()+1, strain_rate_.end()-1, strain_rate_tilde_square);
    auto i = std::distance(strain_rate_.begin()+1, upper);
    visc_ = (viscosity_[i+1] - viscosity_[i])*(strain_rate_tilde_square - strain_rate_[i])/(strain_rate_[i+1] - strain_rate_[i]) + viscosity_[i];
    tau_ = visc_*DXQY::c2Inv + 0.5;


    
    std::valarray<lbBase_t> omega(DXQY::nQ);
    auto tau_inv = 1.0/tau_;
    for (int q=0; q < DXQY::nQ; ++q)
    {
        omega[q] = -tau_inv*(f[q] - feq[q]); 
    }
    
    
    return omega;
}

#endif

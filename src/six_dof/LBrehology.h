#ifndef LBREHOLOGY_H
#define LBREHOLOGY_H

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "../lbsolver/LBglobal.h"

//=====================================================================================
//
//                          P O W E R   L A W   R H E O L O G Y
//
//=====================================================================================

template <typename DXQY>
class PowerLawRheology
{
public:
    PowerLawRheology(const std::string file_name, const lbBase_t K);
    /*template <typename T1, typename T2, typename T3>
     std::valarray<lbBase_t> omegaBGK(const T1 &f,
                                     const lbBase_t& rho,
                                     const T2 &u,
                                     const lbBase_t& u_sq,
                                     const std::valarray<lbBase_t> &cu,
                                     const T3 &force,
                                     const lbBase_t &source); */
    template<typename T1, typename T2, typename T3> 
    lbBase_t tau(const T1 &f, const lbBase_t& rho, const T2 &u, const lbBase_t& u_sq, const std::valarray<lbBase_t> &cu, const T3 &F) const;                                 
    inline lbBase_t getTau() const {
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
    lbBase_t scaleX_;
    std::vector<lbBase_t> tabular_strain_rate_;
    std::vector<lbBase_t> tabular_viscosity_;
};

//                                 PowerLawRheology
//----------------------------------------------------------------------------------- PowerLawRheology
template <typename DXQY>
PowerLawRheology<DXQY>::PowerLawRheology(const std::string file_name, const lbBase_t K)
/* Class constructor, fills the share rate vs viscosity lookup-table:

    nu = K(gamma)^(n-1)

   We have tabulated 
   <n, double>
   <number of entries, int> 
   x  tau

   x = (K/c_s^2)^(1/(n-1))*sqrt(2E_ijE_ij)/(2c_s^2\rho) 
 */ 
{
    std::ifstream data_table_file(file_name, std::ios::in);

    double power_law_exponent;
    data_table_file >> power_law_exponent;
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

    scaleX_ = std::pow(K*DXQY::c2Inv, 1.0/(power_law_exponent - 1.0))*0.5*DXQY::c2Inv;
}

//                                 PowerLawRheology
//----------------------------------------------------------------------------------- tau
template <typename DXQY>
template<typename T1, typename T2, typename T3> 
lbBase_t PowerLawRheology<DXQY>::tau(const T1 &f, const lbBase_t& rho, const T2 &u, const lbBase_t& u_sq, const std::valarray<lbBase_t> &cu, const T3 &F) const
{
    // Calculate E_ij = cicjf - cs^2rho\delta_ij - rho u_iu_j + 0.5(u_iF_j + u_jF_i)
    auto Eij = DXQY::qSumCCLowTri(f);
    lbBase_t EijEij = 0;

    int cnt = 0;
    for (int i = 0; i < DXQY::nD; ++i) {
        for (int j = 0; j < i; ++j) {
            Eij[cnt] += 0.5*(u[i]*F[j] + u[j]*F[i]) - rho*u[i]*u[j];   
            EijEij += 2*Eij[cnt]*Eij[cnt];
            cnt += 1;
        }
        Eij[cnt] += u[i]*F[i] -  rho*(DXQY::c2 + u[i]*u[i]);
        EijEij += Eij[cnt]*Eij[cnt];
        cnt += 1;
    }
    
    lbBase_t x = scaleX_*std::sqrt(2.0*EijEij); // /rho;

    auto upper = std::upper_bound(tabular_strain_rate_.begin()+1, tabular_strain_rate_.end()-1, x);

    auto i = std::distance(tabular_strain_rate_.begin()+1, upper);


    //std::cout << x << " " << i << " " << upper -  tabular_strain_rate_.begin()+1<< std::endl;


    return (tabular_viscosity_[i+1] - tabular_viscosity_[i])*(x - tabular_strain_rate_[i])/(tabular_strain_rate_[i+1] - tabular_strain_rate_[i]) + tabular_viscosity_[i];
}

#endif
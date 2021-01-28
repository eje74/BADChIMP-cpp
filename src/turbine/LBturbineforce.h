#ifndef LBTURBINEFORCE_H
#define LBTURBINEFORCE_H

#include <valarray>
#include "../lbsolver/LBglobal.h"

template <typename T>
inline std::valarray<lbBase_t> rotatingForce(const lbBase_t &w, const lbBase_t &rho, const T &pos, const T &vel)
{
    lbBase_t rhow2 = rho*w*2;
    lbBase_t rhoww = rho*w*w;
    std::valarray<lbBase_t> F {                         0.0, 
                                rhow2*vel[2] + rhoww*pos[1],
                               -rhow2*vel[1] + rhoww*pos[2]
                               };
    return F;
}

template<typename U, typename V>
inline std::valarray<lbBase_t> velRot2Inert(const lbBase_t omega_x, const U &pos, const V &vel_r)
{
    std::valarray<lbBase_t> vel_i(3);
    
    vel_i[0] = vel_r[0];
    vel_i[1] = vel_r[1] - omega_x*pos[2];
    vel_i[2] = vel_r[2] + omega_x*pos[1];
    
    return vel_i;
}

template<typename U, typename V>
inline std::valarray<lbBase_t> velInert2Rot(const lbBase_t omega_x, const U &pos, const V &vel_i)
{
    return velRot2Inert(-omega_x, pos, vel_i);
}


#endif

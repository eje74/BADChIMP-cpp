#ifndef LBSNIPPETS_H
#define LBSNIPPETS_H

#include "LBfield.h"
#include "LBglobal.h"

// Set gravity body force
template <typename DXQY>
inline std::valarray<lbBase_t> setForceGravity(const lbBase_t &rho0, const lbBase_t &rho1, const VectorField<DXQY> &vField, const int &nodeNo)
{
    return (rho0/(rho0 + rho1)) * vField(0, nodeNo);
}

// Constant pressure with oil
inline void setConstDensity(lbBase_t &q0, lbBase_t &q1, lbBase_t &rho0, lbBase_t &rho1, const lbBase_t rhoConst)
// Assuming index 0: water and 1:oil
// Pressure is kept constant by adding oil and removing oil/water mixture
{

    lbBase_t rhoTot = rho0 + rho1;
    lbBase_t qTmp = 2*(rhoConst - rhoTot);

    lbBase_t x1 = 1.0  - (rho0/rhoTot)*(qTmp < 0.0);

    q1 = x1 * qTmp;
    q0 = qTmp - q1;

    rho0 += 0.5*q0;
    rho1 += 0.5*q1;
}

inline void setConstSource(lbBase_t &q0, lbBase_t &q1, lbBase_t &rho0, lbBase_t &rho1, const lbBase_t rate)
// Assuming index 0: water and 1:oil
// adding oil/water mixture at constant rate
{
    q0 = rate * rho0/(rho0 + rho1);
    q1 = rate - q0;
    rho0 += 0.5*q0;
    rho1 += 0.5*q1;
}


inline std::valarray<lbBase_t> sphereCenter(const std::vector<int>& rotCenterPos, const lbBase_t r0, const lbBase_t theta0, const lbBase_t angVel, const int t)
  //
{
  std::valarray<lbBase_t> ret(3);
  ret[0] = rotCenterPos[0] + r0*cos(theta0+ angVel*t);
  ret[1] = rotCenterPos[1] + r0*sin(theta0+ angVel*t);
  ret[2] = rotCenterPos[2];

  return ret;
}


template <typename DXQY, typename T1>
  inline std::valarray<lbBase_t> rotatingPointForce(const T1 &f, const std::vector<int>& x, const std::vector<int>& rotCenterPos, const std::valarray<lbBase_t>& sphereCenter, const lbBase_t r0, const lbBase_t theta0, const lbBase_t angVel, const lbBase_t R, const lbBase_t epsilon, const int t)
//
{
  
  std::valarray<lbBase_t> ret(DXQY::nD);
  for(int i = 0; i < DXQY::nD; ++i){
    ret[i] = 0.0;
  }

  if(sphereCenter[0] != rotCenterPos[0] + r0*cos(theta0+ angVel*t))
    std::cout<<"POSITION OF SPHERE DOES NOT MATCH POSITION OF FORCE"<<std::endl;
  
  lbBase_t pi = 3.1415;

  lbBase_t fromSphereCenterSq = 0.0;
  for(int i = 0; i < DXQY::nD; ++i){
    fromSphereCenterSq += (x[i]-sphereCenter[i])*(x[i]-sphereCenter[i]);
  }
  
  if ( fromSphereCenterSq <= (R+epsilon)*(R+epsilon)){
    lbBase_t fromCentOfRotSq = 0.0;
    for(int i = 0; i < DXQY::nD; ++i){
      fromCentOfRotSq += (x[i]-rotCenterPos[i])*(x[i]-rotCenterPos[i]);
    }
    
    lbBase_t fromCentOfRot = sqrt(fromCentOfRotSq);
    
    std::valarray<lbBase_t> velRot(DXQY::nD);
    
    velRot[0] = - angVel*fromCentOfRot*sin(theta0+ angVel*t);
    velRot[1] = angVel*fromCentOfRot*cos(theta0+ angVel*t);

    for(int i = 2; i < DXQY::nD; ++i){
      velRot[i] = 0.0;
    }
    
    std::valarray<lbBase_t> vel(DXQY::nD);
    vel = velRot;
    lbBase_t delta = 1;

    if (fromSphereCenterSq > R*R && epsilon>0){
    
      lbBase_t fromSphereShell = sqrt(fromSphereCenterSq) - R;

      //delta = (1+cos(pi*fromSphereShell/epsilon))/(2*epsilon);
      delta = (1+cos(pi*fromSphereShell/epsilon))*0.5;
    }
    
    ret = 2*(vel*DXQY::qSum(f)- DXQY::qSumC(f))*delta;
  }
  
  return ret;
}


template <typename DXQY, typename T1>
  inline lbBase_t qSrcConstSphere(const T1 &f, const std::vector<int>& x, const std::valarray<lbBase_t>& sphereCenter, const lbBase_t R, const lbBase_t massTarget, const lbBase_t epsilon)
//
{
  lbBase_t ret = 0;

  lbBase_t fromSphereCenterSq = 0.0;
  for(int i = 0; i < DXQY::nD; ++i){
    fromSphereCenterSq += (x[i]-sphereCenter[i])*(x[i]-sphereCenter[i]);
  }
  
  if (fromSphereCenterSq <= R*R){
    ret = 2*(massTarget-DXQY::qSum(f));
  }
  /*
  if ( fromSphereCenterSq <= (R+epsilon)*(R+epsilon)){
    lbBase_t delta = 1;
    if (fromSphereCenterSq > R*R){
      lbBase_t pi = 3.1415;
      lbBase_t fromSphereShell = sqrt(fromSphereCenterSq) - R;
      delta = (1+cos(pi*fromSphereShell/epsilon))*0.5;
    }
    ret = 2*(massTarget-DXQY::qSum(f))*delta;
  }
  */
  
  return ret;
}

template <typename DXQY>
  inline std::valarray<lbBase_t> sphereShellUnitNormal(const std::vector<int>& pos, const std::valarray<lbBase_t>& sphereCenter)
//
{
  std::valarray<lbBase_t> ret(DXQY::nD);
  for(int i = 0; i < DXQY::nD; ++i){
    ret[i] = 0.0;
  }

  lbBase_t fromSphereCenterSq = 0.0;
  for(int i = 0; i < DXQY::nD; ++i){
    fromSphereCenterSq += (pos[i]-sphereCenter[i])*(pos[i]-sphereCenter[i]);
  }
  lbBase_t fromSphereCenter = sqrt(fromSphereCenterSq);

  if (fromSphereCenter>0.0){
    for(int i = 0; i < DXQY::nD; ++i){
      ret[i] = (pos[i]-sphereCenter[i])/fromSphereCenter;
    }
  }
  
  return ret;
}


template <typename DXQY>
  inline std::valarray<lbBase_t> tangentialVector(const std::valarray<lbBase_t>& vec, const std::valarray<lbBase_t> & unitNormal)
//
{
  std::valarray<lbBase_t> ret(DXQY::nD);
  lbBase_t vecDotUnitNormal = 0;
  for(int i = 0; i < DXQY::nD; ++i){
    vecDotUnitNormal += vec[i]*unitNormal[i];
  }
  for(int i = 0; i < DXQY::nD; ++i){
    ret[i] = vec[i] - vecDotUnitNormal*unitNormal[i];
  }
  
  return ret;
}

template <typename DXQY>
  inline std::valarray<lbBase_t> tangentialUnitVector(const std::valarray<lbBase_t>& vec, const std::valarray<lbBase_t> & unitNormal)
//
{
  std::valarray<lbBase_t> ret(DXQY::nD);
  lbBase_t vecDotUnitNormal = 0;
  lbBase_t tangentialNorm = 0;
  
  for(int i = 0; i < DXQY::nD; ++i){
    vecDotUnitNormal += vec[i]*unitNormal[i];
  }
  for(int i = 0; i < DXQY::nD; ++i){
    ret[i] = vec[i] - vecDotUnitNormal*unitNormal[i];
    tangentialNorm += ret[i]*ret[i];
  }
  tangentialNorm = sqrt(tangentialNorm);
  if(tangentialNorm>0.0){
    for(int i = 0; i < DXQY::nD; ++i){
      ret[i]/=tangentialNorm;
    }
  }
  return ret;
}

void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels);

void printAsciiToScreen(int nX, int nY, ScalarField& val, int** lab, double pauseInSeconds);

void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** lab, double pauseInSeconds);


#endif // LBSNIPPETS_H

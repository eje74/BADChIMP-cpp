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


inline std::valarray<lbBase_t> sphereCenter(const std::vector<int>& centerPos, const lbBase_t r0, const lbBase_t theta0, const lbBase_t angVel, const int t)
  //
{
  std::valarray<lbBase_t> ret(3);
  ret[0] = centerPos[0] + r0*cos(theta0+ angVel*t);
  ret[1] = centerPos[1] + r0*sin(theta0+ angVel*t);
  ret[2] = centerPos[2];

  return ret;
}


template <typename DXQY, typename T1>
  inline std::valarray<lbBase_t> rotatingPointForce(const T1 &f, const std::vector<int>& x, const std::vector<int>& centerPos, const std::valarray<lbBase_t>& sphereCenter, const lbBase_t r0, const lbBase_t theta0, const lbBase_t angVel, const lbBase_t R, const lbBase_t epsilon, const int t)
//
{
  
  std::valarray<lbBase_t> ret = {0.0, 0.0, 0.0};

  if(sphereCenter[0] != centerPos[0] + r0*cos(theta0+ angVel*t))
    std::cout<<"POSITION OF SPHERE DOES NOT MATCH POSITION OF FORCE"<<std::endl;
  
  lbBase_t pi = 3.1415;

  lbBase_t fromSphereCenterSq = (x[0]-sphereCenter[0])*(x[0]-sphereCenter[0]) + (x[1]-sphereCenter[1])*(x[1]-sphereCenter[1]) + (x[2]-sphereCenter[2])*(x[2]-sphereCenter[2]);
  
  if ( fromSphereCenterSq <= (R+epsilon)*(R+epsilon)){
  
    lbBase_t fromCentOfRot = sqrt((x[0]-centerPos[0])*(x[0]-centerPos[0])
				  + (x[1]-centerPos[1])*(x[1]-centerPos[1])
				  + (x[2]-centerPos[2])*(x[2]-centerPos[2]));
    
    std::valarray<lbBase_t> velRot(3);
    velRot[0] = - angVel*fromCentOfRot*sin(theta0+ angVel*t);
    velRot[1] = angVel*fromCentOfRot*cos(theta0+ angVel*t);
    velRot[2] = 0.0;

    std::valarray<lbBase_t> vel(3);
    vel = velRot;
    lbBase_t delta = 1;

    if (fromSphereCenterSq > R*R){
    
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

  lbBase_t fromSphereCenterSq = (x[0]-sphereCenter[0])*(x[0]-sphereCenter[0]) + (x[1]-sphereCenter[1])*(x[1]-sphereCenter[1]) + (x[2]-sphereCenter[2])*(x[2]-sphereCenter[2]);
  
  if (fromSphereCenterSq <= R*R){
    ret = 2*(massTarget-DXQY::qSum(f));
  }

  return ret;
}



void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels);

void printAsciiToScreen(int nX, int nY, ScalarField& val, int** lab, double pauseInSeconds);

void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** lab, double pauseInSeconds);


#endif // LBSNIPPETS_H

# -*- coding: UTF-8 -*-
__author__ = 'olau'

import numpy as np
from scipy import optimize

"""
This script generates a table of viscosity values as a function of the 
contraction Etilde_{ij}Etilde_{ij} (EE). This is related to the 
lattice Boltzmann (LB) strainrate gammadot=sqrt(2*EE)/(2*rho*C_2*tau), where rho
is the LB density, C_2 is the lattice constant, and tau is the LB relaxation 
time.
To generate the chosen rheology, one must give the specified model parameters as
input to the mutest-function. This function takes the rheology function as first
input. 

The generated table is saved in prescribed path
"""
                  
def Carreau(x, EE, mu_0, mu_inf, lambda1, n, y0, rho=1.0):
    """
    x is the dyn. visc., EE is the contraction of LB strain rate tensor
    """   
    C_2=1/3.
    tau=x/(rho*C_2)+0.5
    gammadot=np.sqrt(2*EE)/(2*rho*C_2*tau)
    
    rheo_funk=(mu_inf+(mu_0-mu_inf)*(1+(lambda1*gammadot)**y0)**((n-1)/y0)) 
    
    return x-rheo_funk
                                                 

                                                            
def mutest(rheo_root, *argsIn):
    """
    INPUT
    *argsIn are the arguments passed on to the rheo_funk
    """
    lowerlim=1.e-4
    upperlim=1.e10


    return optimize.brentq(rheo_root, lowerlim, upperlim, args=argsIn, maxiter=200)    
                           
                                                 
def write_table(fileName, muList,EEList):                
    with open(fileName, "w") as f:
        f.write("%d\n" % (muList.shape[0]))
        for z in np.arange(muList.shape[0]):
            f.write("%.10e %.10e\n" % (EEList[z], muList[z]))
    print('data written to ' + fileName)

                        

"""
Carreau fluid model
mu_eff = mu_inf+(mu_0-mu_inf)*(1+(lambda1*gammadot)**y0)**((n-1)/y0)
gammadot: (strainrate) sqrt(2*EE)/(2*rho*C_2*tau)

INPUT
EETest: tabulated values for 2*rho*C_2*tau*gammadot (gammadot: Contracted LB strainrate tensor)
mu_0: viscosity at zero shear
mu_inf: viscosity at infinite shear
lambda1: time prameter determening onset of shear thinning
n: shear thinning index
y0: tuning parameter 

RETURNS solution for mu_eff
"""


#Carreau input-------------------
mu_inf=1.666666666666667e-1 #viscosity at infinite shear
mu_0=1e3*mu_inf             #viscosity at zero shear
lambda1=1e7                 #time prameter determening onset of shear thinning
n=0.5                       #shear thinning index
y0=2                        #tuning parameter 
#--------------------------------
#path to, and name of, rheology table being generated   
directory_path_test = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/PythonScripts/test/"
file_name_base_test = "test"
#Setting number of data points for tabulated viscosity
numDataPoints=1000

#Setting upper and lower limit on EE range
lowLimitPower=-20
topLimitPower=0
#initializing EE array evenly distributed numbers on a log scale
testEE1 = np.logspace(lowLimitPower, topLimitPower, numDataPoints)

file_name = directory_path_test + file_name_base_test + '.dat'

testmu = np.zeros(testEE1.shape, dtype=float)

#Generating viscosity tables
for index, EE in enumerate(testEE1):
    testmu[index] = mutest(Carreau,EE, mu_0, mu_inf, lambda1, n, y0)
                         
#Write to files                                                                                                  
write_table(file_name, testmu,testEE1)    
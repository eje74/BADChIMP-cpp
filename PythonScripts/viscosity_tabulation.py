# -*- coding: UTF-8 -*-
__author__ = 'olau'

#import lbstd as lb
#import diffusion as lb_diff
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from mpl_toolkits.mplot3d.axes3d import get_test_data


def Papanastasiou(x, EE, rho, mu_p, tau_yield, m):
    """
    x is relaxation dyn. visc., EE is the contraction of LB strain rate tensor, B is self.rho
    """
    return HerschelBulkley(x, EE, rho, mu_p, tau_yield, 1, m)
    
def HerschelBulkley(x, EE, rho, mu_p, tau_yield, n, m):
    """
    x is relaxation dyn. visc., EE is the contraction of LB strain rate tensor, B is self.rho
    """
    C_2=1/3.
    tau=x/(rho*C_2)+0.5
    gammadot=np.sqrt(2*EE)/(2*rho*C_2*tau)
    
      
    return x-(mu_p*gammadot**(n-1) +tau_yield/(2*gammadot)*(1-np.exp(-m*gammadot)))    

def mutest_Papanastasiou(EETest, mu_p, tau_yield, m):
    lowerlim=1.e-4
    upperlim=1.e10

    if(EETest==0):
        return mu_p+ m*tau_yield                            
    else:
        return optimize.brentq(Papanastasiou, lowerlim, upperlim, args=(EETest, 1.0, mu_p, tau_yield, m), maxiter=200)
        #return mu_p+ m*tau_yield   
        
def mutest_HerschelBulkley(EETest, mu_p, tau_yield, n, m):
    lowerlim=1.e-4
    upperlim=1.e10

    if(EETest==0):
        return mu_p+ m*tau_yield                            
    else:
        return optimize.brentq(HerschelBulkley, lowerlim, upperlim, args=(EETest, 1.0, mu_p, tau_yield, n, m), maxiter=200)
           
                  
def Carreau(x, EE, rho, mu_0, mu_inf, lambda1, n, y0):
    
    
    C_2=1/3.
    tau=x/(rho*C_2)+0.5
    gammadot=np.sqrt(2*EE)/(2*rho*C_2*tau)
    
    return x-(mu_inf+(mu_0-mu_inf)*(1+(lambda1*gammadot)**y0)**((n-1)/y0))   
                                                
        
def mutest_Carreau(EETest, mu_0, mu_inf, lambda1, n, y0):
    """
    Carreau fluid model
    mu_eff = mu_inf+(mu_0-mu_inf)*(1+(lambda1*gammadot)**y0)**((n-1)/y0)
    gammadot: (strainrate) sqrt(2*EE)/(2*rho*C_2*tau)
    
    INPUT
    EETest: tabulated values for 2*rho*C_2*tau*gammadot (Contracted LB strainrate tensor)
    mu_0: viscosity at zero shear
    mu_inf: viscosity at infinite shear
    lambda1: time prameter determening onset of shear thinning
    n: shear thinning index
    y0: tuning parameter 
    
    RETURNS solution for mu_eff
    """
    lowerlim=1.e-4
    upperlim=1.e10

    if(EETest==0):
        return mu_0                            
    else:
        return optimize.brentq(Carreau, lowerlim, upperlim, args=(EETest, 1.0, mu_0, mu_inf, lambda1, n, y0), maxiter=200)

def Carreau_and_degradation(mu_app, EE,Eflow_dir, rho, mu_0, mu_inf, lambda1, lambda2, n, y0, m):
    C_2=1/3.
    tau=x/(rho*C_2)+0.5
    gammadot=np.sqrt(2*EE)/(2*rho*C_2*tau)
    Eflow_dir2=Eflow_dir/(2*rho*C_2*tau)
    
    return mu_app-(mu_inf+(mu_0-mu_inf)*(1+(lambda1*gammadot)**y0)**((n-1)/y0)+ (lambda2*Eflow_dir2)**(m-1))   
            

                    
                                    
def mutest_Carreau_and_degradation(EETest, Eflow_dirTest, mu_0, mu_inf, lambda1, lambda2, n, y0, m):
    """
    Carreau fluid model
    mu_eff = mu_inf+(mu_0-mu_inf)*(1+(lambda1*gammadot)**y0)**((n-1)/y0)
    gammadot: (strainrate) sqrt(2*EE)/(2*rho*C_2*tau)
    +
    Degradation
    
    INPUT
    EETest: tabulated values for 2*rho*C_2*tau*gammadot (Contracted LB strainrate tensor)
    Eflow_dirTest: tabulated values for 2*rho*C_2*tau* u_i u_j E_ij / u**2 (sum over repeating indices)
    mu_0: viscosity at zero shear
    mu_inf: viscosity at infinite shear
    lambda1: time parameter determening onset of shear thinning
    lambda2: time parameter determening onset of shear thickening
    n: shear thinning index
    y0: tuning parameter 
    m: shear thickening index
    
    RETURNS solution for mu_eff
    """
    lowerlim=1.e-4
    upperlim=1.e10

    if(EETest==0 and Eflow_dirTest==0):
        return mu_0                            
    else:
        return optimize.root(Carreau_and_degradation, lowerlim, upperlim, args=(EETest, 1.0, mu_0, mu_inf, lambda1, n, y0), maxiter=200)
 
           
def write_table(fileName, muList,EEList):                
    with open(fileName, "w") as f:
        f.write("%d\n" % (muList.shape[0]))
        for z in np.arange(muList.shape[0]):
            f.write("%.10e %.10e\n" % (EEList[z], muList[z]))
    print 'data written to '+fileName

#______________________________________________________________________
                        
Dhat=40#160
Utarget=2.5e-2#2.5e-2
B=5.
Re=200.  
M=1.e5 #value used in table printed to file
M2=1.e4

n0=0.5
mu_p0=Utarget*Dhat/Re
tau_yield0=B/Re*Utarget**2                                                                                                                                    
m=M*mu_p0/tau_yield0    
m2=M2*mu_p0/tau_yield0                               
#tau_yield0= 2.25e-05#5.625e-06 
#mu_p0=0.006#0.0045
#m=4.e6#1.2e6#12000000.0

mu_p0=1.666666666666667e-1

tau_p0=3.*mu_p0+0.5
print 'tau_p =', tau_p0
print 'tau_yield =', tau_yield0
EEYield=0.5*((tau_yield0/(2*mu_p0))*2*tau_p0/3)**2

#A = Ehat_ij * Ehat_ij

#----------------------------------------------------------
mu_0=1e3*mu_p0
mu_inf=mu_p0
y0=2
n=0.5#0.75
lambda1=1e7

file_name_base2=str(1)
file_name_base="/home/olau/Programs/LB/D2Q9/2016_11_29/tab_visc/tab_visc"+file_name_base2
file_name2 = file_name_base+"Papanastasiou.dat"
file_name3 = file_name_base+"Herschel-Bulkeley.dat"
file_name4 = file_name_base+"PowerLaw.dat"
file_name5 = file_name_base+"Carreau.dat"

file_name_app = np.empty(5, dtype=object)
file_name_app[0] = "Papanastasiou"     
file_name_app[1] = "Herschel-Bulkeley"
file_name_app[2] = "Papanastasiou2"    
file_name_app[3] = "Carreau"
file_name_app[4] = "Power-Law"

file_name = np.empty(file_name_app.shape[0], dtype=object)
for i in range(0,file_name.shape[0]):
    file_name[i] = file_name_base + file_name_app[i] + '.dat'
    

#file_name[0] += "Papanastasiou.dat"     
#file_name[1] += "Herschel-Bulkeley.dat"
#file_name[2] += "Papanastasiou2.dat"    
#file_name[3] += "Carreau.dat"
#file_name[4] += "PowerLaw.dat"
#___________________________________________________________

 
#testStrain=np.logspace(-15,-3,2000)
#testtau=testStrain.copy()
    
sec1Size=np.int(1e2)#np.int(1e4)  
sec2Size=4*sec1Size
sec3Size=sec1Size
secTotSize= sec1Size + sec2Size + sec3Size
    
testEE1=np.logspace(-20,np.log10(EEYield)-2.5,sec1Size+1)
testEE2=np.logspace(np.log10(EEYield)-2.5,np.log10(EEYield)+2.5,sec2Size+1)
testEE3=np.logspace(np.log10(EEYield)+2.5,0,sec3Size+1)




testEE=np.zeros(secTotSize+1)
testEE[0]=0
testEE[1:sec1Size+1]=testEE1[0:sec1Size]
testEE[sec1Size+1:sec1Size+sec2Size+1]=testEE2[0:sec2Size]
testEE[sec1Size+sec2Size+1:secTotSize+1]=testEE3[0:sec3Size]

#testStrain=np.zeros(10001)
#testStrain=np.logspace(-20,0,10001)
#testStrain[0]=0


testMu = np.ones((50,50,1), dtype=float)
#testtau=testStrain.copy()
testmu=np.zeros(testEE.shape[0])
testmu2=np.zeros(testEE.shape[0])
testmu3=np.zeros(testEE.shape[0])
testmu4=np.zeros(testEE.shape[0])
testmu5=np.zeros(testEE.shape[0])



for index in np.arange(testmu.shape[0]):
    testmu[index] = mutest_Papanastasiou(testEE[index],mu_p0,tau_yield0,m)
    testmu2[index] = mutest_HerschelBulkley(testEE[index],mu_p0,tau_yield0,n0,m)
    testmu3[index] = mutest_Papanastasiou(testEE[index],mu_p0,tau_yield0,m2)
    testmu4[index] = mutest_Carreau(testEE[index],mu_0,mu_inf,lambda1,n,y0)
    testmu5[index] = mutest_HerschelBulkley(testEE[index],mu_p0,0,n0,m) #Power-Law


lambda2= 1e7
m=1.5#2.5          



X=testEE
Y=testEE
X, Y = np.meshgrid(X, Y)  
                   
Z=(mu_inf+(mu_0-mu_inf)*(1+(lambda1*X)**y0)**((n-1)/y0))
Z+= (lambda2*Y)**(m-1)
#Z+= mu_0-mu_0*(1+(lambda2*Y))**(m-1)


#Z[Y>1e-2]=mu_inf

Z2=np.zeros(X.shape[0])
for i in np.arange(X.shape[0]):
    Z2[i] = Z[i,i]

                         
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25) 
#X, Y = np.meshgrid(X, Y)                                                                    
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)                                                
                                    
write_table(file_name[0], testmu,testEE)    
write_table(file_name3, testmu2,testEE)  
write_table(file_name[3], testmu4,testEE)  
write_table(file_name4, testmu5,testEE) 


plt.close('all')
matplotlib.rcParams['text.usetex']=True


fig = plt.figure(1,figsize=plt.figaspect(0.33333))
# set up the axes for the first plot
#ax = fig.add_subplot(1, 3, 1, projection='3d')
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=matplotlib.cm.coolwarm,
#                       linewidth=0, antialiased=False)
#fig.colorbar(surf, shrink=0.5, aspect=10)

ax = fig.add_subplot(1, 3, 1)
ax.loglog(X[0,:],Z[0,:])
ax.title.set_text(r'$u_i u_j E_{ij}/u^2=0$')

ax = fig.add_subplot(1, 3, 2)
ax.loglog(Y[:,0],Z[:,0])
ax.title.set_text(r'$\dot{\gamma}=\sqrt{2E_{ij}E_{ij}}=0$')

ax = fig.add_subplot(1, 3, 3)
ax.loglog(X[0,:],Z2)
ax.title.set_text(r'$u_i u_j E_{ij}/u^2=\sqrt{2E_{ij}E_{ij}}$')

#fig = plt.figure(11,figsize=plt.figaspect(1))
## set up the axes for the first plot
#ax = fig.add_subplot(1, 1, 1, projection='3d')
##ax.plot_wireframe(np.log10(X[1:,1:]), np.log10(Y[1:,1:]), np.log10(Z[1:,1:]), rstride=10, cstride=10)
#surf = ax.plot_surface(np.log10(X[1:,1:]), np.log10(Y[1:,1:]), np.log10(Z[1:,1:]), rstride=1, cstride=1, cmap=matplotlib.cm.hot_r,
#                       linewidth=0, antialiased=False)

fig=plt.figure(2)
plt.pcolormesh(np.log10(X[1:,1:]),np.log10(Y[1:,1:]),np.log10(Z[1:,1:]))
plt.colorbar()

fig=plt.figure(21)
plt.loglog(X[0,1:], Z2[1:])



fig=plt.figure(3)
plt.clf()
ax1 = fig.add_subplot(111)
plt.title('testmu, EEYield = '+np.str(EEYield))
plt.loglog(testEE, testmu)
plt.loglog(testEE, testmu2)
plt.loglog(testEE, testmu3)
plt.loglog(testEE, testmu4)
plt.loglog(testEE, testmu5)
plt.legend(['M='+str(M),'Herschel-Bulkley','M='+str(M2), 'Carreau','Power law'])
ax1.set_xlabel('EE')
ax1.set_ylabel('mu')

fig=plt.figure(4)
plt.clf()
ax1 = fig.add_subplot(111)
plt.title('testmu, EEYield = '+np.str(EEYield))
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu+0.5)/3), testmu)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu2+0.5)/3), testmu2)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu3+0.5)/3), testmu3)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu4+0.5)/3), testmu4)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu5+0.5)/3), testmu5)
plt.legend(['M='+str(M),'Herschel-Bulkley','M='+str(M2), 'Carreau','Power law'])
ax1.set_xlabel('gammadot')
ax1.set_ylabel('mu')

stress1=2*testmu*np.sqrt(2*testEE)/(2*(3*testmu+0.5)/3)
stress2=2*testmu2*np.sqrt(2*testEE)/(2*(3*testmu2+0.5)/3)
stress3=2*testmu3*np.sqrt(2*testEE)/(2*(3*testmu3+0.5)/3)
stress4=2*testmu4*np.sqrt(2*testEE)/(2*(3*testmu4+0.5)/3)
stress5=2*testmu5*np.sqrt(2*testEE)/(2*(3*testmu5+0.5)/3)

fig=plt.figure(5)
plt.clf()
ax1 = fig.add_subplot(111)
#plt.plot(np.sqrt(2*testEE)/(2*(3*testmu+0.5)/3), stress1)
plt.plot(np.sqrt(2*testEE)/(2*(3*testmu2+0.5)/3), stress2)
#plt.plot(np.sqrt(2*testEE)/(2*(3*testmu3+0.5)/3), stress3)
#plt.plot(np.sqrt(2*testEE)/(2*(3*testmu4+0.5)/3), stress4)
plt.plot(np.sqrt(2*testEE)/(2*(3*testmu5+0.5)/3), stress5)
plt.legend(['Herschel-Bulkley','Power law'])
ax1.set_xlabel('gammadot')
ax1.set_ylabel('Shear Stress')

fig=plt.figure(6)
plt.clf()
ax1 = fig.add_subplot(111)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu+0.5)/3), stress1)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu2+0.5)/3), stress2)
#plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu3+0.5)/3), stress3)
#plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu4+0.5)/3), stress4)
plt.loglog(np.sqrt(2*testEE)/(2*(3*testmu5+0.5)/3), stress5)
plt.legend(['M='+str(M),'Herschel-Bulkley','Power law'])
ax1.set_xlabel('gammadot')
ax1.set_ylabel('Shear Stress')

plt.draw()
plt.show()

#___________________________________________________________ 
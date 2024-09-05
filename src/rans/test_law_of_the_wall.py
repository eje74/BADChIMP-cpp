import numpy as np
from matplotlib import pyplot as plt

E = 9.7
k = 0.41
d = 0.5
rho = 1
mu = 1.6e-5
u = 0.185

yp = np.linspace(10, 300, 1000)

A = 1/k
B = E*d/mu


f = lambda x : A*x*np.log(B*x) 
df = lambda x : A*np.log(B*x) + A

u_max = f(mu*300/d)
u_min = f(mu*10/d)

print(u_min, u_max)

plt.figure(0)
us = mu*155/d
plt.plot(0, us, 'k.')
plt.figure(1)
plt.plot(0, f(us) - u, 'k.')



if (u > u_min) and (u < u_max):
    for n in range(1, 10+1):
        dus = -(f(us) - u)/df(us)
        us += dus
        plt.figure(0)
        plt.plot(n, us, 'k.')
        plt.figure(1)
        plt.plot(n, f(us) - u, 'k.')
        
plt.show()
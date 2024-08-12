# Python script documentation and use

## C++ functions defined in main
We have two functions for writing boundary conditions BADCHiMP-input files:  
- ```write_boundary_file```: Write information for the use of a link based boundary condition
- ```write_boundary_file_node```: Write information of the use in a node based boundary condition.

### File format boundary node
```python
"""
    Returns
    -------
    np.array
        One entry for each boundary node
        position               |  x, y    |  0,  1
        distance from boundary |  gamma   |  2
        point0                 | (x0, y0) |  3,  4
        point1                 | (x1, y1) |  5,  6
        point2                 | (x2, y2) |  7,  8
        point3                 | (x3, y3) |  9, 10
        distance to cm interpol|  gamma2  | 11
        weight bnd intersection|  w_a[:4] | 12, 13, 14, 15
        weight interp. points  |  w_b[:4] | 16, 17, 18, 19
        weight boundary node   |  w_c[:4] | 20, 21, 22, 23
        wall normal            | nx, ny   | 24, 25

    Notes
    -----
    Explanatory comments to the return values
    - position: is the x and y coordinate to the boundary node. Is also
       used for identification.
    - distance from boundary: -gamma*normal is the nearest point
       on the surface relative to 'position'.
    - point0-3: gives the coordinates to the four points used for the 
       bi-linear interpolation
    - distance to cm interpol: gamma2*normal is the point nearest to the
       center of the four interpolation points along the normal centered 
       at 'position'.
    - weight bnd intersection: interpolation weights such that
             sum w_a[n]*point_n = (x, y) - gamma*normal
    - weight interp. points: interpolation weights such that
             sum w_b[n]*point_n = (x, y) + gamma2*normal
    - weight boundary node: interpolation weights such that
             sum w_c[n]*point_n = (x, y)
    - wall normal: is the outward pointing wall normal 
"""
``` 

## Algorithms

### Definition of shape functions 
See the [comsol documentation](https://doc.comsol.com/5.3/doc/com.comsol.help.comsol/comsol_api_xmesh.40.4.html) for shape-functions.

### Non-equilibrium functions near the wall  
Close to the wall it is assumed that the stress tensor only has a shear component which implies that the shear rate tensor, $E_{ij}$, only has off diagonal non-zero elements. 
This is related to velocity distributions by
$$ E_{ij} = \frac{1}{2\rho c_s^2\tau}\tilde{E}_{ij}, $$
where
$$ \tilde{E}_{ij} = -\big[\sum_\alpha f^\mathrm{neq}_\alpha Q_{\alpha ij} + \frac{1}{2}\big(u_iF_j + u_jF_i\big)\big], $$
where
$$Q_{\alpha ij} = c_{\alpha i}c_{\alpha j} - c_s^2\delta_{ij}.$$  
It is also assumed that the shear stress close to the wall is constant, hence we need to find a correction, $\delta f_\alpha^w$ to the non-equilibrium velocity distributions so that this will be the case. 
Here we set that
$$  \delta f_\alpha^w = w_\alpha Q_{\alpha ij}\frac{\delta\Pi_{ij}}{2c_s^4}. $$  

We want
$$\tilde{E}_{nn} = \tilde{E}_{tt} = 0 $$
where $n$ and $t$ represents the normal and tangential directions to the wall, respectively.  

We can use the following calculation to obtain  $\delta\Pi_{ij}$
$$ \tilde{E}_{ij} =  -\big[\sum_\alpha (f^\mathrm{neq}_\alpha +  \delta f_\alpha^w)Q_{\alpha ij} + \frac{1}{2}\big(u_iF_j + u_jF_i\big)\big] $$
$$ \tilde{E}_{ij} =  -\big[\sum_\alpha Q_{\alpha ij}f^\mathrm{neq}_\alpha +  \delta\Pi_{ij} + \frac{1}{2}\big(u_iF_j + u_jF_i\big)\big], $$
so that
$$ \delta\Pi_{ij} = -\big[ \tilde{E}_{ij} + \sum_\alpha Q_{\alpha ij}f^\mathrm{neq}_\alpha  + \frac{1}{2}\big(u_iF_j + u_jF_i\big) \big] $$
Further, we have that
$$ Q_{\alpha nn} = Q_{\alpha ij}n_in_j = (c_{\alpha i}n_i)^2 - c_s^2 .$$

#### Non equilibrium wall algorithm  (not stable)
- Calculate $c_{\alpha i}n_i$ and $c_{\alpha i}t_i$, where $t_x = -n_y$ and  $t_y = n_x$.  
- Calculate $\Pi_{nn}$ and $\Pi_{tt}$  
- Calculate $\delta\Pi_{nn}$ and $\delta\Pi_{tt}$ as $\delta\Pi_{nn} = -\big[ \Pi_{nn}  + u^w_nF^w_n\big]$ and $\delta\Pi_{tt} = -\big[ \Pi_{tt}  + u^w_tF^w_t\big]$  
- $\delta\Pi_{ij} = -\big[ \Pi_{nn}  + u^w_nF^w_n\big]n_in_j - \big[ \Pi_{tt}  + t^w_tF^w_n\big]t_it_j$
- $\delta f_\alpha^w = w_\alpha Q_{\alpha ij}\frac{\delta\Pi_{ij}}{2c_s^4} = -\frac{w_\alpha}{2c_s^4}\big(Q_{\alpha nn}(\Pi_{nn}  + u^w_nF^w_n) + Q_{\alpha tt}(\Pi_{tt}  + u^w_tF^w_t)\big)$ 
- $f_\alpha^{w,\mathrm{neq}} = f_\alpha^{\mathrm{neq}} + \delta f_\alpha^w$

####  Non equilibrium wall algorithm  (alternative 2)  
Second order expression for the collision term. 

$\Omega_\alpha \approx w_\alpha\frac{Q_{\alpha ij}}{c_s^2}\rho S_{ij} + \frac{c_{\alpha i}F_i}{c_s^2} $.
From the LB equation we have that 

$\Omega_\alpha = -\frac{1}{\tau}f_\alpha^\mathrm{neq} + w_\alpha\big(1- \frac{1}{2\tau}\big)\frac{c_{\alpha i}F_i}{c_s^2}$ 

$f_\alpha^\mathrm{neq} = -\tau w_\alpha\frac{Q_{\alpha ij}}{c_s^2}\rho S_{ij} - w_\alpha \frac{c_{\alpha i}F_i}{c_s^2}$ 

$S_{ij} = \frac{1}{2}\frac{\partial u_t}{\partial\vec{n}}\big(t_in_j + n_it_j\big)$ 

So that  

$f_\alpha^\mathrm{neq} = -\tau w_\alpha\frac{Q_{\alpha nt}}{c_s^2}\rho \frac{\partial u_t}{\partial\vec{n}} - w_\alpha \frac{c_{\alpha i}F_i}{c_s^2}$




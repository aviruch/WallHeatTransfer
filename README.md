
# External wall heat transfer

## Problem
The wall of a cold room is exposed to solar radiation over a period of time. The wall is made of an insulating material with thermal conductivity of 0.038 W/m-K, specific heat capacity of 70 J/kg-K and density of 120 kg/m<sup>3</sup>. The wall can be considered as infinite in two -dimensions and finite along the thickness of 15 cm. Once the wall is exposed to radiation, the temperature across the wall starts increasing from a uniform initial value of 20˚C across thickness. 

Assume an average constant radiation of 650 W/m<sup>2</sup>; outdoor and indoor air temperature maintained at 27˚C and 12˚C, respectively; average heat transfer coefficient to be 15 W/m<sup>2</sup>-K; formulate the transient heat transfer problem in terms of difference equations with appropriate boundary condition. 


## Dependencies

- Python 
- Numpy 
- Scipy 
- Sympy
- Matplotlib 

```python
A = 1 # cross sectional area of wall element in m^2
L = 0.15  # With of the wall in meter
nx = 6  # number of locations on the wall nodes
dx = L / (nx - 1)  # distance between two consecutive locations
k = 0.30 # thermal conductivity of wall material in W / (m*C)
ro = 120 # Density of material
cp = 70 # specific heat capacity in J / (kg*C)
alpha = k/(ro*cp)  #  thermal diffusivity 
print ("thermal diffusivity",alpha)
q_left = 600 # W/sqm
h = 15 # convective heat transfer coefficient in W / (m^2 * C)
```

Simultaneous equations (at each node)


```python



def solver_implicit(dict):
    eq0 = sym.Eq(q_left + h*(35-t0)+ k*(t1-t0)/dx,k*dx*(t0 - dict[0])/(alpha*2*dt))
    eq1 =sym.Eq(k*(t0-t1)/dx + k*(t2-t1)/dx,k*dx*(t1 - dict[1])/(alpha*2*dt))
    eq2 =sym.Eq(k*(t1-t2)/dx + k*(t3-t2)/dx,k*dx*(t2 - dict[2])/(alpha*2*dt))
    eq3 =sym.Eq(k*(t2-t2)/dx + k*(t4-t3)/dx,k*dx*(t3 - dict[3])/(alpha*2*dt))
    eq4 =sym.Eq(k*(t3-t4)/dx + k*(t5-t4)/dx,k*dx*(t4 - dict[4])/(alpha*2*dt))
    eq5 =sym.Eq(k*(t4-t5)/dx + h*(12-t5),k*dx*(t5 - dict[5])/(alpha*2*dt))
    sol = sym.solve([eq0,eq1,eq2,eq3,eq4,eq5],(t0,t1,t2,t3,t4,t5))
    return sol
```

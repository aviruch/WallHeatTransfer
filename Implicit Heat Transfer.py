import numpy
from matplotlib import pyplot

def solver_implicit(dict):
    eq0 = sym.Eq(q_left + h*(35-t0)+ k*(t1-t0)/dx,k*dx*(t0 - dict[0])/(alpha*2*dt))
    eq1 =sym.Eq(k*(t0-t1)/dx + k*(t2-t1)/dx,k*dx*(t1 - dict[1])/(alpha*2*dt))
    eq2 =sym.Eq(k*(t1-t2)/dx + k*(t3-t2)/dx,k*dx*(t2 - dict[2])/(alpha*2*dt))
    eq3 =sym.Eq(k*(t2-t2)/dx + k*(t4-t3)/dx,k*dx*(t3 - dict[3])/(alpha*2*dt))
    eq4 =sym.Eq(k*(t3-t4)/dx + k*(t5-t4)/dx,k*dx*(t4 - dict[4])/(alpha*2*dt))
    eq5 =sym.Eq(k*(t4-t5)/dx + h*(12-t5),k*dx*(t5 - dict[5])/(alpha*2*dt))


    sol = sym.solve([eq0,eq1,eq2,eq3,eq4,eq5],(t0,t1,t2,t3,t4,t5))


    return sol

def flatSol(sol_t):
    l= []
    l.append(sol_t[t0])
    l.append(sol_t[t1])
    l.append(sol_t[t2])
    l.append(sol_t[t3])
    l.append(sol_t[t4])
    l.append(sol_t[t5])
    return l

# Main Starts here

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
# Define the locations along the rod.
x = numpy.linspace(0.0, L, num=nx)

# Set the initial temperature along the rod.
T0 = numpy.full((10,1),25)
print ("To",T0)
##T0[0] = 35.0 # Outside Temperature K
##T0[-1] = 12.0 # Inside Temperature K

# Set the time-step size based on CFL limit.
nt = 25  # number of time steps to compute
sigma = 0.5
dt = 20 #sigma * dx**2 / alpha  # time-step size

print ("dt",dt)

# Compute the temperature along the wall.

#Node 0

import sympy as sym
#sym.init_printing()
t0,t1,t2,t3,t4,t5 = sym.symbols('t0,t1,t2,t3,t4,t5')

sol_t_20 = solver_explicit([25,25,25,25,25,25,25])
print ("sol_t_20",sol_t_20)

# Solution 40 sec
sol_t_40 = solver_explicit(flatSol(sol_t_20))
print ("sol_t_40",sol_t_40)

# Solution 60 sec
sol_t_60 = solver_explicit(flatSol(sol_t_40))
print ("sol_t_60",sol_t_60)

# Solution 80 sec
sol_t_80 = solver_explicit(flatSol(sol_t_60))
print ("sol_t_80",sol_t_80)

# Solution 100 sec
sol_t_100 = solver_explicit(flatSol(sol_t_80))
print ("sol_t_100",sol_t_100)

# Solution 120 sec
sol_t_120 = solver_explicit(flatSol(sol_t_100))
print ("sol_t_120",sol_t_120)

# Solution 140 sec
sol_t_140 = solver_explicit(flatSol(sol_t_120))
print ("sol_t_140",sol_t_140)

masterlist = []
masterlist.append(flatSol(sol_t_20))
masterlist.append(flatSol(sol_t_40))
masterlist.append(flatSol(sol_t_60))
masterlist.append(flatSol(sol_t_80))
masterlist.append(flatSol(sol_t_100))
masterlist.append(flatSol(sol_t_120))
masterlist.append(flatSol(sol_t_140))
print ("Master",masterlist)

##
##for impacts in masterlist:
##    pyplot.plot(impacts)
##     
##pyplot.show()

#Plot the temperature along the wall.
pyplot.figure(figsize=(6.0, 4.0))
pyplot.xlabel('Node')
pyplot.ylabel('Temperature [C]')
pyplot.grid()
#pyplot.plot(x, flatSol(sol_t_20),"r",flatSol(sol_t_40),"bs", flatSol(sol_t_60),flatSol(sol_t_80))
#pyplot.plot(x, T, color='C0', linestyle='-', linewidth=2)
i=0
for impacts in masterlist:
    i=i+1
    pyplot.plot(impacts,label="{} data".format(i))
    pyplot.legend()

#pyplot.xlim(0.0, L)
#pyplot.ylim(0.0, 100.0);
pyplot.show()









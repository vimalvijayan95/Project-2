# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 17:39:21 2019

@author: andrea
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
#Program to solve by the finite element method de diffusion (heat equation).

def analytic(x,t):
    term1 = np.exp(-(np.pi**2)*t)
    term2 = np.sin(np.pi*x)
    return term1*term2

#---------------------|
# 1st: \Delta x=1/10; |
#-------------------- |   

# Space interval:
L = 1 
xa = 0;
xb = L;
Nx = 100;
x = np.zeros(Nx+1)
x[0] = xa 

# Time interval:
T = 1;
ta = 0;
tb = T;
Nt = 60000
t = np.zeros(Nt+1)
t[0] = ta

#Differentials:
dx = (xb-xa)/Nx;
dt = (tb-ta)/Nt;

# Time and space values:
u = np.zeros((Nx+1,Nt+1))
for i in range(Nx):
    x[i+1] = xa + dx*(i+1)
for j in range(Nt):
    t[j+1] = ta + dt*(j+1)

# Amplitude:
u[:,0] = np.sin(np.pi * x)
# Boundary conditions: 
u[0,:] = 0;
u[Nx,:] = 0;

# Use Euler's method to solve the problem:
#------------------------------------------------
d = dt/(dx**2)

# To make our equation stable we need that s <= 1/2

for j in range(Nt):
    for i in range(1,Nx):
        u[i,j+1] = u[i,j] + d*(u[i+1,j] + u[i-1,j]- 2*u[i,j])
       
# Transpose our matrix to match values with analytical solution:
ut = np.zeros((Nt+1,Nx+1))
for j in range(Nt+1):
    for i in range(Nx+1):
#Compare with the analytical solution:
        ut[j,i] = u[i,j]

# Evaluate the difference between the numerical and analytical solution: 
X,T = np.meshgrid(x,t)
u_ann = analytic(X,T)
diff = np.abs(u_ann - ut)
print('Max. difference between Euler method and analytical solution:',np.max(diff))


# Save the results:
np.save('Euler_600', ut)
np.save('An_Euler_600',u_ann)
#%%

#--------------------------
#EULER METHOD OPTIMIZATION|
#--------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

def analytic(x,t):
    term1 = np.exp(-(np.pi**2)*t)
    term2 = np.sin(np.pi*x)
    return term1*term2

L = 1; 
xa = 0;
xb = L;
Nx = 100;
x = np.zeros(Nx+1)
x[0] = xa

T = 1;
ta = 0;
tb = T;

# Define different values for the time interval. 
Ntv = np.arange(20000,210000,10000)

k = 0;
maxdiff = np.zeros(len(Ntv))

for Nt in Ntv:
    t = np.zeros(Nt+1)
    t[0] = ta
# The difference:
    dx = (xb-xa)/Nx;
    dt = (tb-ta)/Nt;

# Create the different elements we will use:
    u = np.zeros((Nx+1,Nt+1))

    for i in range(Nx):
        x[i+1] = xa+dx*(i+1)
    for j in range(Nt):
        t[j+1] = ta+dt*(j+1)

    u[:,0] = np.sin(np.pi * x)
    u[0,:] = 0;
    u[Nx,:] = 0;
    
    d = dt/(dx**2)

    for j in range(Nt):
        for i in range(1,Nx):
            u[i,j+1] = u[i,j] + d*(u[i+1,j] + u[i-1,j] - 2*u[i,j])
      
    ut = np.zeros((Nt+1,Nx+1))
    for j in range(Nt+1):
        for i in range(Nx+1):
            ut[j,i] = u[i,j]

    X,T = np.meshgrid(x,t)
    u_ann = analytic(X,T)
    diff = np.abs(u_ann - ut)
    maxdiff[k] = np.max(diff)
    print(k,np.max(diff))
    k+=1

np.save('maxdiff_100', maxdiff)
np.save('timesteps_100', Ntv)
#indxm=np.where(maxdiff == np.min(maxdiff))
#3print(Ntv[indxm])


#!/usr/bin/env python
# coding: utf-8

# ## 12 Steps to Navier-Stokes Equations

# ## STEP 1 :       1-D Linear Convection

# Lets start by importing necessary libraries 

# In[94]:


import numpy as np
import matplotlib.pyplot as plt
from math import *


# Defining the domain for both Space & Time 

# In[95]:


xmax = 2                         # Total Domain lenth in (m)
nx   = 41                       # Number of Grid Points in space 'X'
dx   = xmax/(nx-1)               # Size of each grid in (m)
nt   = 25                        # Number of Grid Points in time 't'
dt   = 0.025                     # Time-step size in (s)
c    = 1                         # Wave Propagation Speed in (m/s)
x    = np.linspace(0,xmax,nx)    # Space Domain with grids in (m)


# In[96]:


print(x)


# Defining the Initial Conditions 

# In[97]:


u = np.ones(nx)                      # Initialise the velocity array of ones
u[int(.5/dx):int((1/dx)+1)] = 2      # Implementing the square wave condition for velocity
print(u)


# Plotting the Velocity field function U to understand its profile variation in Space

# In[98]:


plt.plot(x,u, label = 'Initial Velocity Profile')
plt.title('Velocity Profile - Square Wave')
plt.xlabel('X (m)')
plt.ylabel('u (m/s)')
plt.grid()


# Discretisation of Linear Convection equation & finding its max grid resolution to avoid solution blow-up  
# 
# $$\frac{\partial u}{\partial t} + c* \frac{\partial u}{\partial x} = 0$$ 
# 
# Forward Differencing in Time :
# $$\frac{\partial u}{\partial t} \approx \frac{u_i^{n+1}-u_i^{n}}{\Delta t}\rightarrow 1$$
# 
# Backward Differencing in Space :
# $$\frac{\partial u}{\partial x} \approx \frac{u_i^{n}-u_{i-1}^{n}}{\Delta x}\rightarrow 2$$
# 
# Therefore expressing Velocity explicitly, we have;
# $$u_i^{n+1} = u_i^{n} - c*\frac{\Delta t}{\Delta x}*(u_i^{n}-u_{i-1}^{n}) \rightarrow 3$$  
# 
# where $c*\frac{\Delta t}{\Delta x}$ is called the Stability Criterion of the FDM Scheme.
# 
# The CFL Criteria for the given Schemes
# 
# CFL Number :
# $$ c*\frac{\Delta t}{\Delta x} \leq 1 $$
# 
# Hence the conditions for the max grid resolutions are below:
# nx = 81;
# dx = 0.025;
# 
# hence CFL = 1
# 

# Lets start by Intialising a new array of velocity $u_n$

# In[99]:


#Initialize a temporary array
un = np.ones(nx) 

# Time Loop
for n in range(nt):
    un = u.copy()
# Space Loop
    for i in range(1,nx):
        u[i] = un[i]-c*dt/dx*(un[i]-un[i-1])


# Lets Plot the following 

# In[100]:


plt.plot(x,u,label='Convected Velocity Profile')
plt.title('Convected Profile Velocity')
plt.xlabel('X (m)')
plt.ylabel('u (m/s)')
plt.grid()
plt.show


# In[ ]:





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 13:00:08 2021

@author: louis
"""

import numpy as np
import scipy as sp
import networkx as nx
import scipy.integrate as spode
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

import graphSimple as gs
import graphProcess as gp


k_Al = 167
cp_Al = 896 #J/(kg.C-1)
rho_Al = 2700 # Kg/m**3

L_beam = 0.3 #m
A_beam = 0.01 #m^2
M_beam = rho_Al*L_beam*A_beam

R_beam = L_beam/(k_Al*A_beam)
C_beam = cp_Al*M_beam


################ Using a classical thermal evolution scheme


# n_lumps = 10

# C_lump = C_beam/n_lumps
# R_lump = R_beam/n_lumps


# simG = gp.SimGraph()
# gs.defineGraphBeamSim(n_lumps, conductance = 1/R_lump, Capa = C_lump, SimGraph = simG)
# simG.modifyNodeType(n_lumps-1, 'const')
# simG.plotGraph()
# simG.prepareSimulation()

# thermal_ode = lambda T, time: simG.thermalConductionDirect(T, time)

# # Profile parabolique
# x = np.linspace(0, 1, n_lumps)
# T_max = 500
# T_min = 300
# T0 = T_min + (T_min-T_max)*4*x*(x-1)
# t = np.linspace(0, 100, 100)

# test_evo = simG.thermalConductionDirect(T0, t)
# r=np.multiply(test_evo,simG.adj_mat)


# sol = spode.odeint(thermal_ode, T0, t)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X, Y = np.meshgrid(np.arange(n_lumps), t)
# Z = sol

# ax.plot_surface(X, Y, Z,cmap=cm.coolwarm, linewidth=0, antialiased=False)

# ax.set_xlabel('Lump Number')
# ax.set_ylabel('Time')
# ax.set_zlabel('Temperature')

# plt.show()


############### Using the pre-computed matrix 


n_lumps = 10

C_lump = C_beam/n_lumps
R_lump = R_beam/n_lumps


simG = gp.SimGraph()
gs.defineGraphBeamSim(n_lumps, conductance = 1/R_lump, Capa = C_lump, SimGraph = simG)
simG.modifyNodeType(n_lumps-1, 'const')
simG.plotGraph()
simG.prepareSimulation()

thermal_ode = lambda T, time: simG.thermalConductionLin(T, time)

# Profile parabolique
x = np.linspace(0, 1, n_lumps)
T_max = 500
T_min = 300
T0 = T_min + (T_min-T_max)*4*x*(x-1)

lt = np.linspace(0, 100, 100)

sol = spode.odeint(thermal_ode, T0, lt)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(np.arange(n_lumps), lt)
Z = sol

ax.plot_surface(X, Y, Z,cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_xlabel('Lump Number')
ax.set_ylabel('Time')
ax.set_zlabel('Temperature')

plt.show()


simG.addNode('var', C_lump)
simG.addEdge(simG.n_lumps - 2, simG.n_lumps - 1, edge_value= 2, edge_type=1)
simG.plotGraph()

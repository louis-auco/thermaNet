#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 21:38:29 2021

@author: louis
"""

import numpy as np
import networkx as nx

# Defining the data needed in a graph for the simulation
# Each node will either be a variable or a constant, with actions defined on the nodes
# At this point no physical domain change is to be made, with the first versions 
# destined to do thermal simulations.
# It is meant to simulate thermal conduction in lumped modelsdefined through graphs.
# The constant nodes will be used as a mean to simulate Dirichlet conditions, with 
# null Neumann conditions considered as the default on all nodes.
# The Neumann conditions will be added through an additional if not null.

class SimGraph:
    # "Defining a graph for simulation"
    def __init__(self):
        self.graph = nx.Graph()
        self.n_lumps = 0
        

    def addNode(self, sim_type, capacitance = 0):
        if sim_type in ['const', 'var']:
            self.graph.add_node(self.n_lumps, sim_type = sim_type, capacitance = capacitance, actions = [])
            self.n_lumps += 1
            
        else:
            print('sim_type must be equal to const or var')
    
    def modifyNodeType(self, node_id, sim_type):
        if sim_type in ['const', 'var']:
            self.graph.nodes[node_id]['sim_type']=sim_type
        else:
            print('sim_type must be equal to const or var')
    
    def addEdge(self, node0, node1, edge_value, edge_type = 0):
        # In the thermal conduction case, the edge value is the conductance between the lumps
        # Consider edge_type 0 to be conduction and 1 radiation
        try:
            self.graph.add_edge(node0, node1, edge_value = edge_value, edge_type = edge_type)
        except:
            print("Can't add this edge (already existing)")
    
    def modifyEdgeType(self, edge_tuple, edge_type):
        node0 = edge_tuple[0]
        node1 = edge_tuple[1]
        self.graph[node0][node1]['edge_type'] = edge_type
    
    def modifyEdgeValue(self, edge_tuple, edge_value):
        node0 = edge_tuple[0]
        node1 = edge_tuple[1]
        self.graph[node0][node1]['edge_value'] = edge_value
        
    def addAction(self, list_nodes, action):
        if hasattr(list_nodes, '__iter__'):
            node_set = list(set(list_nodes))
            for i in node_set:
                self.graph.nodes[i]['actions'].append(action)
        else:
            self.graph.nodes[list_nodes]['actions'].append(action)
    
    
    def plotGraph(self):
        edges = self.graph.edges()
        def catchColor(edge):
            try:
                return edge['color']
            except:
                return 'b'
        
        cmap = { 0:'b',1:'y',2:'k',3:'g',4:'r' }
        edgeColors = [cmap[self.graph[u][v]['edge_type']] for u,v in edges]
        nodeColors = ['g' if node[1]['sim_type']=='var' else 'r' for node in self.graph.nodes(data=True)]
        nx.draw(self.graph, with_labels=True, edge_color=edgeColors,node_color=nodeColors, font_weight='bold')
    
    
    def MatAdj_type(self):
        self.adj_type = nx.adjacency_matrix(self.graph, dtype='float32',weight='edge_type').todense()
        return self.adj_type
    
    
    def Mat_Radiation(self):
        self.adj_type = nx.adjacency_matrix(self.graph, dtype='float32',weight='edge_type').todense()
        self.adj_rad = np.multiply(self.adj_mat, self.adj_type == 1)
        return self.adj_rad
     
       
    def Mat_Conduction(self):
        self.adj_mat = nx.adjacency_matrix(self.graph, dtype='float32',weight='edge_value').todense()
        self.adj_cond = np.multiply(self.adj_mat, self.adj_type == 0)
        return self.adj_cond
    
    
    def eigenSimplify(self, mat_problem):
        ### Return the eigenvalues of the thermal conduction problem
        n_tot = self.n_lumps
        # Eigenvalue computations
        # The eigenvalue computation cannot be made directly through the Mat_R
        # matrix, derived from the adjacency matrix. The following computation
        # allows to reconstruct the matrix
        res = np.zeros((n_tot,n_tot,n_tot))
        e = np.eye(n_tot)
        
        O1 = np.ones((n_tot,1))
        OT = O1.T
        
        # Computing a basis for the 1(n,n) - 1(n,n)' application
        # res is a n_tot*n_tot*n_tot matrix
        for i in range(n_tot):
            for j in range(n_tot):
                temp_e = np.reshape(e[:,i], (self.n_lumps,1))
                temp = O1 @ temp_e.T - temp_e @ OT
                res[i,:,:] = temp
        
        # Multiply each basis matrix of the res matrix by the resistance matrix
        r_mat = np.zeros((n_tot,n_tot,n_tot))
        for i in range(n_tot):
            r_mat[i,:,:] = np.multiply( np.squeeze(res[i,:,:]), mat_problem)

        # Summing the matrices on the 3rd dimension
        mat_t = np.zeros(n_tot)
        for i in range(n_tot):
            mat_t = mat_t + np.squeeze(r_mat[:,:,i])

        C_m1 = np.reshape(np.power(self.Capa_vect,-1), (self.n_lumps,1))
        # Getting the real matrix by incorporating the thermal capacitances
        mat_t = (C_m1@O1.T)*mat_t
        # Getting the eigenvalues
        eig_R = np.linalg.eigvals(mat_t)
        return eig_R, mat_t
        # mat_t is used after for the computation of the thermal evolution
        
        
    def defineCapacitance(self):
        self.Capa_vect = np.array(list(nx.get_node_attributes(self.graph,'capacitance').values()))
        self.sim_type_vect = list(nx.get_node_attributes(self.graph,'sim_type').values())
        self.const_nodes = [i for i,sim_type in enumerate(self.sim_type_vect) if sim_type == 'const']
        for i in self.const_nodes:
            self.Capa_vect[i]=np.inf


    def prepareSimulation(self):
        self.MatAdj_type()
        self.Mat_Conduction()
        self.defineCapacitance()
        self.eig_Cond, self.mat_Cond = self.eigenSimplify(self.adj_cond)
        self.eig_Rad, self.mat_Rad = self.eigenSimplify(self.adj_rad)
        
    def thermalConductionDirect(self, T, time):
        
        # Using the adjacency matrix the computations of the fluxes can be
        # simplified
        # mat_T is a matrix for which each column contains the same temperature Ti
        mat_T = np.tile(T, (self.n_lumps,1))
        # mat_dT is an antisymmetric matrix containing all possible temperature
        # differences Ti-Tj
        mat_dT = mat_T.T-mat_T
        # Using the Hadamard product, between mat_dT and the resistance matrix
        # computed thanks to the adjacency matrix, we get the thermal fluxes
        # between the lumps
        Q_cond_mat = np.multiply(mat_dT,self.adj_mat)
        # We sum each row to get the total flux for each lump and then apply
        # the capacitance vector and the external power (B_vect) to get the final dT
        Q_cond_vect = np.ravel(np.sum(Q_cond_mat,0))
        dT = np.multiply(np.power(self.Capa_vect,-1),(Q_cond_vect))
        return dT
    
    def thermalRadiationLin(self, T, time):
        # Using the simplified matrix to compute the thermal evolution caused by radiation
        dT = self.mat_Rad @ np.power(T, 4)
        return dT
    
    def thermalConductionLin(self, T, time):
        # Using the simplified matrix to compute the thermal evolution caused by diffusion
        dT = self.mat_Cond @ T
        return dT
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 13:00:36 2021

@author: louis
"""
import networkx as nx

def plotGraph(graph):
    edges = graph.edges()
    colors = [graph[u][v]['color'] for u,v in edges]
    nx.draw(graph, with_labels=True, edge_color=colors, font_weight='bold')
    
    
def defineGraphBeam(n_lumps, conductance, Capa):
    graph = nx.Graph()
    graph.add_node(0, Capa = Capa)
    for i in range(n_lumps-1):
        graph.add_node(i+1, Capa = Capa)
        graph.add_edge(i,i+1, color='b', conductance = conductance)
    return graph

    
def defineGraphBeamSim(n_lumps, conductance, Capa, SimGraph):
    SimGraph.addNode(sim_type = 'const', capacitance = Capa)
    for i in range(n_lumps-1):
        SimGraph.addNode(sim_type = 'var', capacitance = Capa)
        SimGraph.addEdge(i,i+1, edge_value = conductance)

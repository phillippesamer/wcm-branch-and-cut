#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy
import networkx as nx
import matplotlib.pyplot as plt

import os

# five G(n,p) random graphs, for p in {0.01, 0.02, ..., 0.25}
n = 50
percentual_rates = range(1, 26, 1)      # q \in 1, 2, ..., 25 => p \in 0.01, 0.02, ... , 0.25
num_examples = 5
make_figures = False
destination_folder = "./new_g_np"

if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

seeds = [12182243, 45266141]
random.seed(seeds[0])
numpy.random.seed(seeds[1])

for q in percentual_rates:
    
    for iteration in range(num_examples):
        
        p = q / 100
        G = nx.erdos_renyi_graph(n, p)
        m = G.number_of_edges()
        
        vertex_weights = random.sample(range(-50, 51), n)
        
        filepath = f"{destination_folder}/{n}-{q}-{iteration}"
        
        # save the graph to png file
        if make_figures:
            fig = plt.figure()
            nx.draw(G, with_labels=True, node_color=[x+50 for x in vertex_weights], node_size=300, cmap=plt.cm.Blues)
            fig.savefig( filepath + ".png" )
            plt.close(fig)
        
        # create the .gcc file
        file = open(filepath + ".gcc", "w")
        file.write(f"# G_{n,p} (Erdos-Renyi) graph, with n = {n} and p = {p}\n")
        file.write(f"# Example {iteration} of {num_examples}\n")
        
        file.write(f"{filepath}\n")
        file.write(f"{n}\n")
        file.write(f"{m}\n")

        for i,j in G.edges:
            file.write(f"{i} {j}\n")

        for u in range(n):
            file.write(f"{vertex_weights[u]}\n")

        file.close()

"""
# generate random number of colours as well?
import random
random.seed(621077887)
num_components = [random.randint(3,14) for i in range(25)]
for x in num_components:
    for i in range(5):
        print(x)
"""

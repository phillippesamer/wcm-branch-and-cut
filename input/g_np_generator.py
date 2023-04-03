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
destination_folder = "./alternate_g_np"

if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

seeds = [1092593, 1984337]
random.seed(seeds[0])
numpy.random.seed(seeds[1])

for q in percentual_rates:
    
    for iteration in range(num_examples):
        
        p = q / 100
        G = nx.erdos_renyi_graph(n, p)
        m = G.number_of_edges()
        
        #vertex_weights = random.sample(range(1, 101), n)        # all different
        vertex_weights = [random.randint(1,100) for u in range(n)]
        
        filepath = f"{destination_folder}/{n}-{q}-{iteration}"
        
        # save the graph to png file
        if make_figures:
            fig = plt.figure()
            nx.draw(G, with_labels=True, node_color=[x+50 for x in vertex_weights], node_size=300, cmap=plt.cm.Blues)
            fig.savefig( filepath + ".png" )
            plt.close(fig)
        
        # create the .gcc file
        file = open(filepath + ".gcc", "w")
        file.write(f"# G_{n,p} (Erdos-Renyi) graph, with n = {n} and p = {p}, and non-negative vertex weights\n")
        file.write(f"# Example {iteration} of {num_examples}\n")
        
        file.write(f"{filepath}\n")
        file.write(f"{n}\n")
        file.write(f"{m}\n")

        for i,j in G.edges:
            file.write(f"{i} {j}\n")

        for u in range(n):
            file.write(f"{vertex_weights[u]}\n")

        file.close()

# number of colours (constant for each combination of n,p)
num_rates = len(percentual_rates)
num_components = [random.randint(2,11) for i in range(num_rates)]

file = open(destination_folder + "_colours.txt", "w")

for x in num_components:
    for i in range(5):
        file.write(f"{x}\n")

file.close()

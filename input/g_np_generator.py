#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy
import networkx as nx
import matplotlib.pyplot as plt

import os

# five G(n,p) random graphs, for p in {0.01, 0.02, ..., 0.10}
n = 1000
percentual_rates = range(1, 11, 1)      # q \in 1, 2, ..., 10 => p \in 0.01, 0.02, ... , 0.10
num_examples = 5
make_figures = False
destination_folder = "./g_np_1000"
min_weight = 1
max_weight = 100

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
        
        edge_weights = [random.randint(min_weight,max_weight) for e in range(m)]
        # all different (NB! set min/max weights above accordingly)
        #edge_weights = random.sample(range(min_weight,max_weight), m)
        
        filepath = f"{destination_folder}/{n}-{q}-{iteration}"
        
        # save the graph to png file
        if make_figures:
            fig = plt.figure()
            nx.draw(G,
                    with_labels = True,
                    edge_color = [x for x in edge_weights],
                    node_size = 300,
                    cmap = plt.cm.Blues)
            fig.savefig( filepath + ".png" )
            plt.close(fig)
        
        # create the .gcc file
        file = open(filepath + ".gcc", "w")
        file.write(f"# G_{n,p} (Erdos-Renyi) graph, with n = {n} and p = {p}, and non-negative edge weights\n")
        file.write(f"# Example {iteration} of {num_examples}\n")
        
        file.write(f"{filepath}\n")
        file.write(f"{n}\n")
        file.write(f"{m}\n")

        idx = 0
        for i,j in G.edges:
            file.write(f"{i} {j} {edge_weights[idx]}\n")
            idx = idx+1

        file.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy
import networkx as nx
import matplotlib.pyplot as plt

import os

# five G(n,p) random graphs, for p in {0.01, 0.02, ..., 0.10}
n = 5000
#percentual_rates = range(1, 11, 1)      # q \in 1, 2, ..., 10 => p \in 0.01, 0.02, ... , 0.10
percentual_rates = range(5, 51, 5)      # q \in 5, 10, ..., 50 => p \in 0.05, 0.1, ... , 0.5
num_examples = 5
make_figures = False
destination_folder = f"./g_np_{n}"

gaussian_weights = True  # gaussian or uniform random weights?
weight_mean = 50.0
weight_stddev = 50.0
min_weight = -1000       # only used if gaussian_weights
max_weight = 1000        # is set to False

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
        
        if gaussian_weights:
            edge_weights = [round(random.gauss(mu=weight_mean, sigma=weight_stddev)) for e in range(m)]
            min_weight = min(edge_weights)
            max_weight = max(edge_weights)
        else:
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
        file.write(f"# G_{n,p} (Erdos-Renyi) graph, with n = {n} and p = {p}, and edge weights in [{min_weight}, {max_weight}]\n")
        file.write(f"# Example {iteration} of {num_examples}\n")
        
        file.write(f"{filepath}\n")
        file.write(f"{n}\n")
        file.write(f"{m}\n")

        idx = 0
        for i,j in G.edges:
            file.write(f"{i} {j} {edge_weights[idx]}\n")
            idx = idx+1

        file.close()

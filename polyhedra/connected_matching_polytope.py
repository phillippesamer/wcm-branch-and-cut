#!/usr/bin/env python
# coding: utf-8

from separators import min_ab_separators

import networkx as nx
import subprocess

def get_graph_P5() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,6))
    G.add_edges_from([(1,2), (2,3), (3,4), (4,5)])
    return G

def get_graph_C4() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,5))
    G.add_edges_from([(1,2), (2,3), (3,4), (4,1)])
    return G

def get_graph_claw() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,5))
    G.add_edges_from([(1,2), (1,3), (1,4)])
    return G

def get_graph_paw() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,5))
    G.add_edges_from([(1,2), (1,3), (2,3), (3,4)])
    return G

def get_graph_diamond() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,5))
    G.add_edges_from([(1,2), (1,3), (1,4), (2,3), (3,4)])
    return G

def get_graph_k_1_4() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,6))
    G.add_edges_from([(1,2), (1,3), (1,4), (1,5)])
    return G

def get_graph_k_1_5() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,7))
    G.add_edges_from([(1,2), (1,3), (1,4), (1,5), (1,6)])
    return G

def get_ilp_formulation(G: nx.Graph) -> str:
    """
    Create a string representing the .lp file for the 'separators formulation'
    of connected matchings in G
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()

    # cutset of each vertex
    cutset = [[] for u in range(n+1)]
    for idx, (v1, v2) in enumerate(G.edges()):
        cutset[v1].append(idx+1)
        cutset[v2].append(idx+1)

    # strings corresponding to variable names
    vars = {}
    # TO DO: iterate edges of the graph?
    for e in range(1,m+1):
        vars[e] = " x" + str(e) + " "

    # arbitrary objective function, with constant coefficients
    objective = "Maximize\n  "
    for e in range(1,m+1):
        objective += vars[e]
        if e < m:
            objective += "+ "
        else:
            objective += "\n"

    constraints = "Subject To\n"

    # 1. degree inequalities
    for u in range(1,n+1):
        constraints += "DEG" + str(u) + ": "
        for count, edge_idx in enumerate(cutset[u]):
            constraints += "x" + str(edge_idx)
            if count < len(cutset[u])-1:
                constraints += " + "
            else:
                constraints += " "
        constraints += "<= 1\n"

    # 2. separator inequalities
    for u in range(1,n+1):
        for v in range(u+1, n+1):
            if (G.has_edge(u,v) == False):
                print("sep("+str(u)+","+str(v)+")")
                separators = min_ab_separators(G,u,v)
                for idx,Z in enumerate(separators):
                    #label
                    constraints += "(" + str(u) + "," + str(v) + ")-SEP_#" + str(idx+1) + ": "

                    #x_u + x_v - \sum_(z in Z) x_z <= 1
                    for count, edge_idx in enumerate(cutset[u]):
                        constraints += "x" + str(edge_idx)
                        if count < len(cutset[u])-1:
                            constraints += " + "
                        else:
                            constraints += " "

                    for edge_idx in cutset[v]:
                        constraints += "+x" + str(edge_idx)

                    for z in Z:
                        for count, edge_idx in enumerate(cutset[z]):
                            constraints += " -x" + str(edge_idx)

                    constraints += " <= 1\n"

    # not applicable in this formulation
    bounds = "Bounds\n"

    domain = "Binaries\n"
    for e in range(1,m+1):
        domain += vars[e]
    domain += "\n"
 
    sections = [objective, constraints, bounds, domain, "End\n"]
    contents = "\n".join(sections)
    return contents

def main():
    # specify the input graph following the example functions above
    graph = get_graph_claw()
    output_lp_file     = "examples/claw_original.lp"
    output_facets_file = "examples/claw_facets.lp"

    # generate "separators-based formulation" ILP corresponding to this input
    print("writing ilp on file " + output_lp_file)

    ilp = get_ilp_formulation(graph)

    f = open(output_lp_file, 'w')
    f.write(ilp)
    f.close()

    print("done\n")

    print("NB! For ease of reference, the edge list is")
    for idx, (v1, v2) in enumerate(graph.edges()):
        print(f"edge #{idx+1} is ({v1}, {v2})")

    # run polymake to enumerate facets of the convex hull of integer points above
    print("\nrunning polymake and writing facets on file " + output_facets_file)

    proc = subprocess.run(
    ["polymake", "--script", "polymake_cvxhull_from_lp.pl", output_lp_file],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True)

    f = open(output_facets_file, 'w')
    f.write(proc.stdout)
    f.close()

    print(proc.stderr)
    print("done\n")


if __name__ == "__main__":
    main()

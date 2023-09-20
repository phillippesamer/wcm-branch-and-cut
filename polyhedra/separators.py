#!/usr/bin/env python3
# coding: utf-8

import networkx as nx
import matplotlib.pyplot as plt

# forked from https://github.com/v-kam/shen_liangli_1997
def min_ab_separators(G: nx.Graph, a, b):
    """
    Implementation of algorithm described in "Efficient enumeration of all
    minimal separators in a graph" by H. Shen, W. Liangi on Theoretical
    Computer Science 180 (1997) pp. 169-180

    :param G: input graph
    :param a: start node
    :param b: end node
    :return: all minimal ab-separators (list of sets) 
    """
    n = len(G)
    if type(a) is not set:
        a = {a}
    if type(b) is not set:
        b = {b}

    def N(X):
        """Return set of neighbours of vertices in subset X"""
        if type(X) is not set:
            X = {X}
        neighbors = set()
        for v in X:
            neighbors.update(nx.all_neighbors(G, v))
        return neighbors - X

    def I(X, Cb):
        if type(X) is not set:
            X = {X}
        isolated = set()
        for v in X:
            if v not in N(Cb):
                isolated.add(v)
        return isolated

    def calculate_Cb(S):
        H = G.copy()
        H.remove_nodes_from(S)

        for C in nx.connected_components(H):
            if len(b - C) == 0:
                return C
        return set()

    Cb = calculate_Cb(N(a))
    k = 0
    L = {k: [(a, N(a) - I(N(a), Cb))]}
    for i in range(1, n - 2):
        L[i] = []
    separators = [N(a) - I(N(a), Cb)]

    while k <= n - 3 and len(Cb) != 0:
        for p, S in L[k]:
            for x in S:
                if len(b.intersection(N(x))) == 0:
                    N_plus = N(x) - a - S
                    SuN_plus = S.union(N_plus)
                    Cb = calculate_Cb(SuN_plus)

                    if len(Cb) != 0:
                        S_new = SuN_plus - I(SuN_plus, Cb)
                        if S_new not in separators:
                            L[k + 1].append((x, S_new))
                            separators.append(S_new)
        k += 1
#        if not L[k]:
        if k not in L:
            break

    if separators == [set()]:
        return []
    else:
        return separators

def test_min_ab_separators_paper_example() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,8))

    G.add_edges_from([(1,3), (1,4), (3,4)])
    G.add_edges_from([(1,5), (3,7), (4,6)])
    G.add_edges_from([(2,6),(2,7),(5,6), (5,7)])

    sep = min_ab_separators(G, 1, 2)

    # set of sets to make the comparison easy, ignoring order of separators
    assert {frozenset(s) for s in sep} == \
    { frozenset({5, 3, 6}), \
      frozenset({4, 5, 7}), \
      frozenset({6, 7}), \
      frozenset({5, 4, 3}) }

def test_min_ab_separators_tcsstack_example() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,8))

    G.add_edges_from([(1,2), (1,3)])
    G.add_edges_from([(2,3)])
    G.add_edges_from([(2,4), (2,5), (2,6)])
    G.add_edges_from([(3,4), (3,5), (3,6)])
    G.add_edges_from([(4,7), (5,7), (6,7)])

    sep = min_ab_separators(G, 1, 7)

    # set of sets to make the comparison easy, ignoring order of separators
    assert {frozenset(s) for s in sep} == \
    { frozenset({2, 3}), \
      frozenset({4, 5, 6}) }

def test_min_ab_separators_grid() -> nx.Graph:
    G = nx.grid_graph(dim=(2,3))
    sep = min_ab_separators(G, (0,0), (2,1))

    # set of sets to make the comparison easy, ignoring order of separators
    assert {frozenset(s) for s in sep} == \
    { frozenset({ (0,1),(1,0) }), \
      frozenset({ (1,0),(1,1) }), \
      frozenset({ (1,1),(2,0) }) }

def test_min_ab_separators_P5() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1,6))
    G.add_edges_from([(1,2), (2,3), (3,4), (4,5)])

    sep_1_2 = min_ab_separators(G, 1, 2)
    sep_1_3 = min_ab_separators(G, 1, 3)
    sep_1_4 = min_ab_separators(G, 1, 4)
    sep_1_5 = min_ab_separators(G, 1, 5)
    sep_2_3 = min_ab_separators(G, 2, 3)
    sep_2_4 = min_ab_separators(G, 2, 4)
    sep_2_5 = min_ab_separators(G, 2, 5)
    sep_3_4 = min_ab_separators(G, 3, 4)
    sep_3_5 = min_ab_separators(G, 3, 5)
    sep_4_5 = min_ab_separators(G, 4, 5)

    # set of sets to make the comparison easy, ignoring order of separators
    assert {frozenset(s) for s in sep_1_2} == set()
    assert {frozenset(s) for s in sep_1_3} == { frozenset({2}) }
    assert {frozenset(s) for s in sep_1_4} == { frozenset({2}), frozenset({3}) }
    assert {frozenset(s) for s in sep_1_5} == \
      { frozenset({2}), frozenset({3}), frozenset({4}) }
    assert {frozenset(s) for s in sep_2_3} == set()
    assert {frozenset(s) for s in sep_2_4} == { frozenset({3}) }
    assert {frozenset(s) for s in sep_2_5} == { frozenset({3}), frozenset({4}) }
    assert {frozenset(s) for s in sep_3_4} == set()
    assert {frozenset(s) for s in sep_3_5} == { frozenset({4}) }
    assert {frozenset(s) for s in sep_4_5} == set()

def main():
    test_min_ab_separators_paper_example()
    test_min_ab_separators_tcsstack_example()
    test_min_ab_separators_grid()
    test_min_ab_separators_P5()

    #nx.draw(G, with_labels=True)
    #plt.show()
    #plt.savefig("graph.png")

if __name__ == "__main__":
    main()

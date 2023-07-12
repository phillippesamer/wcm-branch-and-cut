#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <iostream>
#include <limits>

#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/core.h>

using namespace std;
using namespace lemon;

/***
 * \file graph.h
 * 
 * Module containing different data structures representing the same graph:
 * - an edge list, with terminals of each edge stored in vectors s and t
 * - an adjacency matrix storing the corresponding index (or -1 for non-edges)
 * - an adjacency list from the LEMON (Library for Efficient Modeling and 
 * Optimization in Networks), so as to use the highly efficient implementations
 * of algorithms they offer
 * 
 * Some classes are declared friends to avoid cumbersome get/set calls.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 02.07.2023
 */
class Graph
{
public:
    Graph();
    Graph(long, long);
    virtual ~Graph();

    void init_index_matrix();
    void free_index_matrix();

    void init_lemon();

private:
    friend class IO;
    friend class WCMModel;
    friend class WCMCutGenerator;

    long num_vertices;
    long num_edges;

    vector<double> w;    // edge weights

    // adjacency list
    vector< list<long> > adj_list;

    // edge list
    vector<long> s;      // terminal vertex 1
    vector<long> t;      // terminal vertex 2

    // adjacency matrix storing edge indexes
    bool using_matrix;
    long **index_matrix;

    // LEMON adjacency list: http://lemon.cs.elte.hu/pub/doc/1.3/a00237.html
    bool using_lemon;
    ListGraph *lemon_graph;
    vector<ListGraph::Node> lemon_vertices;
    vector<ListGraph::Edge> lemon_edges;
    ListGraph::EdgeMap<long> *lemon_edges_inverted_index;
    ListGraph::EdgeMap<double> *lemon_weight;

    bool lemon_test_adj(ListGraph &, ListGraph::Node &, ListGraph::Node &);
    ListGraph::Edge lemon_test_adj_getting_edge(ListGraph &,
                                                ListGraph::Node &,
                                                ListGraph::Node &);
};

#endif

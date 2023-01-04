#include "graph.h"

Graph::Graph()
{
    using_matrix = false;
    using_lemon = false;
    num_vertices = 0;
    num_edges = 0;
}

Graph::Graph(long n, long m)
{
    using_matrix = false;
    using_lemon = false;
    num_vertices = n;
    num_edges = m;

    adj_list.reserve(n);
    adj_list.insert( adj_list.begin(), n, list<long>() );

    s.reserve(m);
    t.reserve(m);

    w.reserve(n);
}

Graph::~Graph()
{
    adj_list.clear();

    s.clear();
    t.clear();
 
    w.clear();

	if (using_matrix)
        free_index_matrix();

    if (using_lemon)
    {
        lemon_vertices.clear();
        lemon_edges.clear();
        delete lemon_weight;
        delete lemon_edges_inverted_index;
        delete lemon_graph;
    }
}

void Graph::init_index_matrix()
{
	using_matrix = true;

    index_matrix = new long*[num_vertices];
    for (long i=0; i<num_vertices; ++i)
    {
        index_matrix[i] = new long[num_vertices];
        for (long j=0; j<num_vertices; ++j)
            index_matrix[i][j] = -1;
    }
}

void Graph::free_index_matrix()
{
    for (long i = 0; i<num_vertices; ++i)
        delete[] index_matrix[i];

    delete[] index_matrix;
}

void Graph::init_lemon()
{
    using_lemon = true;

    lemon_graph = new ListGraph();

    lemon_vertices.reserve(num_vertices);
    for (long i=0; i<num_vertices; ++i)
        lemon_vertices.push_back(lemon_graph->addNode());

    lemon_edges.reserve(num_edges);
    lemon_edges_inverted_index = new ListGraph::EdgeMap<long>(*lemon_graph);

    // populate lemon graph from the edge list in this object
    for (long idx=0; idx<num_edges; ++idx)
    {
        long v1 = s.at(idx);
        long v2 = t.at(idx);

        ListGraph::Edge e = lemon_graph->addEdge(lemon_vertices[v1], lemon_vertices[v2]);
        lemon_edges.push_back(e);
        (*lemon_edges_inverted_index)[e] = idx;
    }

    // vertex to weight map
    lemon_weight = new ListGraph::NodeMap<double>(*lemon_graph);
    for (long idx=0; idx<num_vertices; ++idx)
    {
        ListGraph::Node v = lemon_vertices[idx];
        double weight = w.at(idx);
 
        (*lemon_weight)[v] = weight;
    }
}

ListGraph::Edge Graph::lemon_test_adj_getting_edge(ListGraph &g,
                                                   ListGraph::Node &x,
                                                   ListGraph::Node &y)
{
    /***
     * Auxiliary function to test adjacency in the LEMON data structure.
     * Returns edge, if found; otherwise, returns INVALID
     */

    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return e;

    return INVALID;
}

bool Graph::lemon_test_adj(ListGraph &g,
                           ListGraph::Node &x,
                           ListGraph::Node &y)
{
    /// auxiliary function to test adjacency in the LEMON data structure

    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return true;

    return false;
}

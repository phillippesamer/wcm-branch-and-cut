# Connected matching polytope facets 
Enumeration of facets of the connected matching (CM) polytope of a graph: the convex hull of characteristic vectors of edge subsets inducing a matching, whose covered vertices induce a connected subgraph.

### Dependencies

We use [polymake](https://polymake.org) and python>=3.8, with the networkx library.

### How to use this

To inspect facets of the CM polytope of _G_, one should

- add a function on `connected_matching_polytope.py` to create the desired input graph _G_ (just mimic one of the examples in the beginning, _e.g._ `get_graph_Petersen()`)
- edit the first lines of the `main()` function to set the desired output file paths and use the graph of choice, _e.g._ graph = get_graph_Petersen()
- run `python connected_matching_polytope.py`

Have fun! (:


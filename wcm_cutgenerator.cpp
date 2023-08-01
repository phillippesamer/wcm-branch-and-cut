#include "wcm_cutgenerator.h"

/// algorithm setup switches

bool SEPARATE_MSI = true;              // MSI = minimal separator inequalities
bool SEPARATE_BLOSSOM = true;
bool SEPARATE_INDEGREE = true;

bool MSI_STRATEGY_FIRST_CUT_BELOW_ROOT = true;
bool MSI_FROM_INTEGER_POINTS_ONLY = false;

bool BLOSSOM_AT_ROOT_ONLY = false;
bool BLOSSOM_HEURISTIC_SEPARATION = true;

bool INDEGREE_AT_ROOT_ONLY = true;
bool MSI_ONLY_IF_NO_INDEGREE = false;

// clean any bits beyond the corresponding precision to avoid numerical errors?
// (at most 14, since gurobi does not support long double yet...)
// NB! THIS OPTION MIGHT RISK MISSING A VIOLATED INEQUALITY
const bool CLEAN_VARS_BEYOND_PRECISION = false;
const int SEPARATION_PRECISION = 14;

// use a tolerance in dealing with relaxation values?
// e.g. y_u < epsilon instead of y_u == 0
// set to 0 to stick to exact values and operators
const double MSI_EPSILON = 1e-5;
const double MSI_ZERO = MSI_EPSILON;
const double MSI_ONE = 1.0 - MSI_EPSILON;
const double INDEGREE_EPSILON = 1e-5;

///////////////////////////////////////////////////////////////////////////////

/// specialized depth-first search to identify/count connected components

void inline dfs_to_tag_components(long u,
                                  long count, 
                                  vector<long> &components, 
                                  vector< vector<long> > &adj_list)
{
    // auxiliary dfs to check connected components in adj_list

    components[u] = count;

    const long degree = adj_list[u].size();

    for (long i=0; i < degree; ++i)
    {
        long v = adj_list[u].at(i);
        if (components[v] < 0)
            dfs_to_tag_components(v, count, components, adj_list);
    }
}

long inline check_components(vector< vector<long> > &adj_list,
                             vector<long> &components)
{
    /// NB! Expects components initialized as vector<long>(adj_list.size(), -1)

    long count = 0;
    const long n = adj_list.size();

    for (long u=0; u < n; ++u)
    {
        if (components[u] < 0)
        {
            dfs_to_tag_components(u, count, components, adj_list);
            ++count;
        }
    }

    return count;
}

///////////////////////////////////////////////////////////////////////////////

bool inline check_integrality(double *point, long dim)
{
    /// check if all coordinates of a point are zero/one-valued

    for (long i=0; i < dim; ++i)
    {
        if (point[i] > MSI_ZERO && point[i] < MSI_ONE)
            return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////

WCMCutGenerator::WCMCutGenerator(GRBModel *model, GRBVar *x_vars, GRBVar *y_vars, IO *instance)
{
    this->model = model;
    this->x_vars = x_vars;
    this->y_vars = y_vars;
    this->x_integral = false;
    this->y_integral = false;
    this->instance = instance;
    this->num_vertices = instance->graph->num_vertices;
    this->num_edges = instance->graph->num_edges;

    this->at_root_relaxation = true;

    this->blossom_counter = 0;
    this->indegree_counter = 0;
    this->minimal_separators_counter = 0;
    this->msi_next_source = 0;

    /***
     * Support graph (using LEMON) to separate blossom inequalities (BI)
     * We construct the support graph only once, and update only the edge
     * capacities from the current relaxation values.
     */

    this->bi_support_graph = new ListGraph();
    this->bi_support_vertices.reserve(this->num_vertices + 1);
    this->bi_support_edges.reserve(this->num_edges + this->num_vertices);
    this->bi_support_capacity = NULL;

    if (SEPARATE_BLOSSOM)
    {
        // n vertices from instance graph, plus an artificial/dummy universal one
        for (long i=0; i < num_vertices+1; ++i)
            bi_support_vertices.push_back(bi_support_graph->addNode());

        // m edges as in the instance graph, plus one from dummy to each vertex
        for (long idx=0; idx<num_edges; ++idx)
        {
            long v1 = instance->graph->s.at(idx);
            long v2 = instance->graph->t.at(idx);
            ListGraph::Edge e = bi_support_graph->addEdge(bi_support_vertices.at(v1),
                                                          bi_support_vertices.at(v2));
            bi_support_edges.push_back(e);
        }

        for (long idx=0; idx<num_vertices; ++idx)
        {
            long v1 = num_vertices;  // dummy vertex
            long v2 = idx;
            ListGraph::Edge e = bi_support_graph->addEdge(bi_support_vertices.at(v1),
                                                          bi_support_vertices.at(v2));
            bi_support_edges.push_back(e);
        }
    }
}

WCMCutGenerator::~WCMCutGenerator()
{
    this->bi_support_vertices.clear();
    this->bi_support_edges.clear();
    delete bi_support_graph;
}

void WCMCutGenerator::callback()
{
    /***
     * The actual callback method within the solver, searching for and adding 
     * dynamically the most violated (if any) indegree inequality, minimal 
     * separator inequality (MSI), and/or blossom inequality.
     */

    try
    {
        // callback from the search at a given MIP node - may include USER CUTS
        if (where == GRB_CB_MIPNODE)
        {
            // node relaxation solution must be available at the current node
            if (this->getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL)
                return;

            // flag when done with the root node relaxation
            if (this->at_root_relaxation)   // initially true
                if (getDoubleInfo(GRB_CB_MIPNODE_NODCNT) > 0)
                    this->at_root_relaxation = false;

            // retrieve relaxation solution
            x_val = this->getNodeRel(x_vars, num_edges);
            y_val = this->getNodeRel(y_vars, num_vertices);
            x_integral = check_integrality(x_val, num_edges);
            y_integral = check_integrality(y_val, num_vertices);

            if (SEPARATE_BLOSSOM && !x_integral)
            {
                if (at_root_relaxation || !BLOSSOM_AT_ROOT_ONLY)
                    run_blossom_separation(ADD_USER_CUTS);
            }

            bool separated = false;
            if (SEPARATE_INDEGREE)
            {
                if (at_root_relaxation || !INDEGREE_AT_ROOT_ONLY)
                    separated = run_indegree_separation(ADD_USER_CUTS);
            }

            if (SEPARATE_MSI)
            {
                if (y_integral || !MSI_FROM_INTEGER_POINTS_ONLY)
                {
                    if (!MSI_ONLY_IF_NO_INDEGREE || !separated)
                    {
                        if (CLEAN_VARS_BEYOND_PRECISION)
                            clean_vars_beyond_precision(SEPARATION_PRECISION);

                        run_minimal_separators_separation(ADD_LAZY_CNTRS);
                    }
                }
            }

            delete[] x_val;
            delete[] y_val;
        }

        // callback from a new MIP incumbent: only LAZY CONSTRAINTS
        else if (where == GRB_CB_MIPSOL)
        {
            // retrieve solution
            y_val = this->getSolution(y_vars, num_vertices);
            y_integral = check_integrality(y_val, num_vertices);

            if (SEPARATE_MSI)
            {
                if (CLEAN_VARS_BEYOND_PRECISION)
                    clean_vars_beyond_precision(SEPARATION_PRECISION);

                run_minimal_separators_separation(ADD_LAZY_CNTRS);
            }

            delete[] y_val;
        }
    }
    catch (GRBException e)
    {
        cout << "Error " << e.getErrorCode()
             << " during WCMCutGenerator::callback(): ";
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Unexpected error during WCMCutGenerator::callback()" << endl;
    }
}

bool WCMCutGenerator::separate_lpr()
{
    /// Interface to be used when solving the LP relaxation only.

    // avoid cut policy distinguishing root node relaxation from others
    this->at_root_relaxation = false;

    try
    {
        // retrieve relaxation solution
        x_val = new double[num_edges];
        for (long e = 0; e < num_edges; ++e)
            x_val[e] = x_vars[e].get(GRB_DoubleAttr_X);

        y_val = new double[num_vertices];
        for (long u = 0; u < num_vertices; ++u)
            y_val[u] = y_vars[u].get(GRB_DoubleAttr_X);

        x_integral = check_integrality(x_val, num_edges);
        y_integral = check_integrality(y_val, num_vertices);

        bool blossom_cut = false;
        bool indegree_cut = false;
        bool msi_cut = false;

        if (SEPARATE_BLOSSOM && !x_integral)
            blossom_cut = run_blossom_separation(ADD_STD_CNTRS);

        if (SEPARATE_INDEGREE && !INDEGREE_AT_ROOT_ONLY)
            indegree_cut = run_indegree_separation(ADD_STD_CNTRS);

        if (SEPARATE_MSI)
        {
            if (!MSI_ONLY_IF_NO_INDEGREE || !indegree_cut)
            {
                if (CLEAN_VARS_BEYOND_PRECISION)
                    clean_vars_beyond_precision(SEPARATION_PRECISION);

                msi_cut = run_minimal_separators_separation(ADD_STD_CNTRS);
            }
        }

        bool separated = (blossom_cut || indegree_cut || msi_cut);

        // last resort: attempt exact BI separation (and indegree..)
        if (!separated && !x_integral && BLOSSOM_HEURISTIC_SEPARATION)
        {
            BLOSSOM_HEURISTIC_SEPARATION = false;
            separated = run_blossom_separation(ADD_STD_CNTRS);
            BLOSSOM_HEURISTIC_SEPARATION = true;

            if (!separated && !y_integral)
                separated = run_indegree_separation(ADD_STD_CNTRS);
        }

        // clean up
        delete[] x_val;
        delete[] y_val;

        return separated;
    }
    catch (GRBException e)
    {
        cout << "Error " << e.getErrorCode() << " during separate_lpr(): ";
        cout << e.getMessage() << endl;
        return false;
    }
    catch (...)
    {
        cout << "Unexpected error during separate_lpr()" << endl;
        return false;
    }
}

void WCMCutGenerator::clean_vars_beyond_precision(int precision)
{
    /// prevent floating point errors by ignoring digits beyond given precision
    for (long u = 0; u < num_vertices; ++u)
    {
        double tmp = y_val[u] * std::pow(10, precision);
        tmp = std::round(tmp);
        y_val[u] = tmp * std::pow(10, -precision);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool WCMCutGenerator::run_blossom_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    if (BLOSSOM_HEURISTIC_SEPARATION)
        // heuristic separation = attempt to find a violated BI quickly, but might fail - runtime in O(n + m)
        model_updated = separate_blossom_heuristically(cuts_lhs, cuts_rhs);
    else
        // exact separation = either find a violated BI, or decide that none exists - runtime in O(n^3 \sqrt(m))
        model_updated = separate_blossom_exactly(cuts_lhs, cuts_rhs);

    if (model_updated)
    {
        // add cuts
        for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
        {
            ++blossom_counter;

            if (kind_of_cut == ADD_USER_CUTS)
                addCut(cuts_lhs[idx] <= cuts_rhs[idx]);

            else if (kind_of_cut == ADD_LAZY_CNTRS)
                addLazy(cuts_lhs[idx] <= cuts_rhs[idx]);

            else // kind_of_cut == ADD_STD_CNTRS
                model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);
        }
    }

    return model_updated;
}

bool WCMCutGenerator::separate_blossom_exactly(vector<GRBLinExpr> &cuts_lhs,
                                               vector<long> &cuts_rhs)
{
    /***
     * Run classical separation algorithm from Padberg & Rao (1982) - still the
     * most efficient for the uncapacitated case cf. "Odd minimum cut sets and
     * b-matchings revisited", 2008, by [Letchford, Reinelt, Theis].
     */

    // 1. DETERMINE UPDATED EDGE CAPACITIES FROM THE CURRENT RELAXATION

    bi_support_capacity = new ListGraph::EdgeMap<double>(*bi_support_graph);

    // edge uv from the instance graph: capacity[uv] = x*_uv
    for (long idx=0; idx < num_edges; ++idx)
    {
        ListGraph::Edge edge = bi_support_edges.at(idx);
        (*bi_support_capacity)[edge] = x_val[idx];
    }

    // edge ru from dummy vertex to u (NB! using x~y linking constraints here!):
    // capacity[ru]  =  1 - \sum_{v neighbour of u} x*_uv  =  1 - y*_u
    for (long idx=0; idx < num_vertices; ++idx)
    {
        long dummy_edge_idx = num_edges + idx;
        ListGraph::Edge edge = bi_support_edges.at(dummy_edge_idx);

        // y_u is the sum of x_uv for v neighbours of u
        (*bi_support_capacity)[edge] = 1.0 - y_val[idx];
    }

    // 2. CONSTRUCT GOMORY-HU CUT TREE OF THE SUPPORT GRAPH

    // this is the runtime bottleneck: O(n^3 sqrt(m)) in this implementation
    GomoryHu<ListGraph, ListGraph::EdgeMap<double> > cut_tree(*bi_support_graph,
                                                              *bi_support_capacity);
    cut_tree.run();

    // 3. LOOK FOR VIOLATED BI FROM MIN-CUT VALUES AT EACH EDGE IN THE CUT TREE

    ListGraph::Node dummy = bi_support_vertices.back();

    // 3.1 TRAVERSE CUT TREE EDGES BY QUERYING THE PREDECESSOR OF EACH VERTEX (EXCEPT THE ROOT) 
    for (long idx=0; idx<num_vertices+1; ++idx)
    {
        ListGraph::Node s = bi_support_vertices.at(idx);
        ListGraph::Node t = cut_tree.predNode(s);
        if(t != INVALID)   // not the cut tree root (n+1 vertices => n edges)
        {
            // 3.2 MINCUT INDUCED BY THIS EDGE OF THE CUT TREE MAY GIVE A VIOLATED
            // BI IF ITS VALUE IS < 1 AND ONE OF THE CUTSETS IS OF ODD CARDINALITY 
            if (cut_tree.predValue(s) < MSI_ONE)
            {
                // 3.2 SIDE 1: HANDLE INDUCED BY THE CUTSET CONTAINING S
                bool cutset_with_s = true;
                long cutset_size = 0;
                vector<long> cutset_vertices = vector<long>();
                vector<bool> cutset_mask = vector<bool>(num_vertices, false);

                for(GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutNodeIt it(cut_tree, s, t, cutset_with_s); it != INVALID; ++it)
                {
                    // ignore the dummy vertex
                    long vertex_id = bi_support_graph->id(it);
                    if (vertex_id != bi_support_graph->id(dummy))
                    {
                        ++cutset_size;
                        cutset_vertices.push_back(vertex_id);
                        cutset_mask.at(vertex_id) = true;
                    }
                }

                if (cutset_size % 2 == 1)
                {
                    long bi_rhs = (cutset_size - 1) / 2;

                    // determine edges with both endpoints in the cutset and check for violation
                    GRBLinExpr violated_constr = 0;
                    double current_lhs = bi_lhs_from_handle(cutset_vertices, cutset_mask, violated_constr);

                    if (current_lhs > bi_rhs)
                    {
                        cuts_lhs.push_back(violated_constr);
                        cuts_rhs.push_back(bi_rhs);
                    }

                    #ifdef DEBUG_BI
                        if (current_lhs > bi_rhs)
                            cout << "### ADDED BI: (...) = " << current_lhs << " > " << bi_rhs << endl;
                    #endif
                }

                // 3.2 SIDE 2: REPEAT FOR THE HANDLE INDUCED BY THE CUTSET CONTAINING T
                cutset_with_s = false;
                cutset_size = 0;
                cutset_vertices.clear();
                cutset_mask = vector<bool>(num_vertices, false);
                for(GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutNodeIt it(cut_tree, s, t, cutset_with_s); it != INVALID; ++it)
                {
                    // ignore the dummy vertex
                    long vertex_id = bi_support_graph->id(it);
                    if (vertex_id != bi_support_graph->id(dummy))
                    {
                        ++cutset_size;
                        cutset_vertices.push_back(vertex_id);
                        cutset_mask.at(vertex_id) = true;
                    }
                }

                if (cutset_size % 2 == 1)
                {
                    long bi_rhs = (cutset_size - 1) / 2;

                    // determine edges with both endpoints in the cutset and check for violation
                    GRBLinExpr violated_constr = 0;
                    double current_lhs = bi_lhs_from_handle(cutset_vertices, cutset_mask, violated_constr);

                    if (current_lhs > bi_rhs)
                    {
                        cuts_lhs.push_back(violated_constr);
                        cuts_rhs.push_back(bi_rhs);
                    }

                    #ifdef DEBUG_BI
                        if (current_lhs > bi_rhs)
                            cout << "### ADDED BI: (...) = " << current_lhs << " > " << bi_rhs << endl;
                    #endif
                }
            }
        }
    }

    delete bi_support_capacity;

    return (cuts_lhs.size() > 0);
}

double WCMCutGenerator::bi_lhs_from_handle(vector<long> &handle_vertices,
                                           vector<bool> &handle_mask,
                                           GRBLinExpr &constr)
{
    /***
     * Determines the edges induced by a given handle, and fills the constraint
     * lhs with the corresponding x_vars. Returns the x_val in the current
     * relaxation.
     */

    double current_lhs = 0.0;
    long handle_size = handle_vertices.size();

    // Choose faster option: m tests whether an edge is induced by the handle
    // vs. h(h-1)/2 adjacency tests within handle vertices only (h = handle_size)
    if ( num_edges < handle_size*(handle_size-1)/2 )
    {
        for (long edge_idx = 0; edge_idx < num_edges; ++edge_idx)
        {
            long v1 = instance->graph->s.at(edge_idx);
            if (handle_mask.at(v1))
            {
                long v2 = instance->graph->t.at(edge_idx);
                if (handle_mask.at(v2))
                {
                    constr += (x_vars[edge_idx]);
                    current_lhs += x_val[edge_idx];
                }
            }
        }
    }
    else
    {
        for (long i = 0; i < handle_size; ++i)
        {
            for (long j = i+1; j < handle_size; ++j)
            {
                long v1 = handle_vertices.at(i);
                long v2 = handle_vertices.at(j);
                long edge_idx = instance->graph->index_matrix[v1][v2];
                if (edge_idx >= 0)
                {
                    constr += (x_vars[edge_idx]);
                    current_lhs += x_val[edge_idx];
                }
            }
        }
    }

    return current_lhs;
}

bool WCMCutGenerator::separate_blossom_heuristically(vector<GRBLinExpr> &cuts_lhs,
                                                     vector<long> &cuts_rhs)
{
    /***
     * Simple heuristic to try to find a violated BI in linear time:
     * 1. Let H be the support graph induced only from vars in (0,1)
     * 2. Let H_i for i in [p] denote the connected components of H
     * 3. For i in [p], if |V(H_i)| is odd, inspect the corresponding BI for violation
     */

    // 1. DETERMINE VERTICES AND EDGES INDUCED BY FRACTIONAL x* ONLY
    
    vector<bool> vertex_mask = vector<bool>(num_vertices, false);
    vector<bool> edge_mask = vector<bool>(num_edges, false);
    get_fractional_info(vertex_mask, edge_mask);

    // 2. DFS TO CHECK CONNECTED COMPONENTS INDUCED BY FRACTIONAL-VALUED EDGES

    vector<bool> seen = vector<bool>(num_vertices, false);
    for (long source = 0; source < num_vertices; ++source)
    {
        if (vertex_mask.at(source) && !seen.at(source))
        {
            vector<long> component_vertices = vector<long>();
            vector<bool> component_mask = vector<bool>(num_vertices, false);

            dfs_from_frac_x_only(edge_mask,
                                 seen,
                                 source,   // only arg. passed by value
                                 component_vertices,
                                 component_mask);

            long size = component_vertices.size();
            double bi_rhs = (size - 1) / 2;

            if (size % 2 == 1)
            {
                // 3. CHECK IF CURRENT ODD COMPONENT AS HANDLE GIVES A BI VIOLATED AT x*
                GRBLinExpr violated_constr = 0;

                // besides checking violation, lift inequality to include x_vars for induced edges at 0
                double current_lhs = bi_lhs_from_handle(component_vertices, component_mask, violated_constr);

                if (current_lhs - bi_rhs > MSI_ZERO)
                {
                    cuts_lhs.push_back(violated_constr);
                    cuts_rhs.push_back(bi_rhs);

                    #ifdef DEBUG_BI
                        cout << "### ADDED BI: (...) = " << current_lhs << " > " << bi_rhs << endl;
                    #endif
                }
            }
        }
    }

    return (cuts_lhs.size() > 0);
}

void inline WCMCutGenerator::get_fractional_info(vector<bool> &vertex_mask,
                                                 vector<bool> &edge_mask)
{
    /// determine support graph "mask" induced by fractional x vars only

    for (long idx=0; idx < num_edges; ++idx)
    {
        if (x_val[idx] > MSI_ZERO && x_val[idx] < MSI_ONE)
        {
            edge_mask.at(idx) = true;

            long v1 = instance->graph->s.at(idx);
            long v2 = instance->graph->t.at(idx);

            vertex_mask.at(v1) = true;
            vertex_mask.at(v2) = true;
        }
    }
}

void inline WCMCutGenerator::dfs_from_frac_x_only(vector<bool> &edge_mask,
                                                  vector<bool> &seen,
                                                  long u,   // NB! only arg. passed by value
                                                  vector<long> &component_vertices,
                                                  vector<bool> &component_mask)
{
    /// auxiliary dfs over subgraph induced by fractional x vars
    seen.at(u) = true;
    component_vertices.push_back(u);
    component_mask.at(u) = true;

    for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
        it != instance->graph->adj_list.at(u).end(); ++it)
    {
        long v = *it;
        long edge_idx = instance->graph->index_matrix[u][v];

        if (edge_mask.at(edge_idx) && !seen.at(v))
        {
            dfs_from_frac_x_only(edge_mask,
                                 seen,
                                 v,   // NB! only arg. passed by value
                                 component_vertices,
                                 component_mask);
         }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool WCMCutGenerator::run_indegree_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    /* run separation algorithm from "On imposing connectivity constraints in
     * integer programs", 2017, by [Wang, Buchanan, Butenko]
     */
    model_updated = separate_indegree(cuts_lhs, cuts_rhs);

    if (model_updated)
    {
        // add cuts
        for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
        {
            ++indegree_counter;

            if (kind_of_cut == ADD_USER_CUTS)
                addCut(cuts_lhs[idx] <= cuts_rhs[idx]);

            else if (kind_of_cut == ADD_LAZY_CNTRS)
                addLazy(cuts_lhs[idx] <= cuts_rhs[idx]);

            else // kind_of_cut == ADD_STD_CNTRS
                model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);
        }
    }

    return model_updated;
}

bool WCMCutGenerator::separate_indegree(vector<GRBLinExpr> &cuts_lhs,
                                        vector<long> &cuts_rhs)
{
    /// Solve the separation problem for indegree inequalities

    vector<long> indegree = vector<long>(num_vertices, 0);

    // 1. COMPUTE INDEGREE ORIENTING EDGES ACCORDING TO RELAXATION SOLUTION
    for (long idx = 0; idx < num_edges; ++idx)
    {
        long u = instance->graph->s.at(idx);
        long v = instance->graph->t.at(idx);
        if (y_val[u] > y_val[v] + INDEGREE_EPSILON)
            indegree.at(v) += 1;
        else
            indegree.at(u) += 1;
    }

    // 2. EVALUATE LHS (1-d[u])*y[u]
    double lhs_sum = 0.;
    for (long u = 0; u < num_vertices; ++u)
        lhs_sum += ( (1 - indegree.at(u)) * y_val[u] );

    // 3. FOUND MOST VIOLATED INDEGREE INEQUALITY (IF ANY) IF LHS > 1
    if (lhs_sum > 1 + INDEGREE_EPSILON)
    {
        // store inequality (caller method adds it to the model)
        GRBLinExpr violated_constr = 0;

        for (long u = 0; u < num_vertices; ++u)
            violated_constr += ( (1 - indegree.at(u)) * y_vars[u] );

        cuts_lhs.push_back(violated_constr);
        cuts_rhs.push_back(1);
    }

    return (cuts_lhs.size() > 0);
}

////////////////////////////////////////////////////////////////////////////////

bool WCMCutGenerator::run_minimal_separators_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    if (y_integral)
    {
        /* run natural separation algorithm based on depth-first search in the
         * support graph - see for example "Thinning out Steiner trees - a node- 
         * based model for uniform edge costs", 2016, by [Fischetti, Leitner,
         * Ljubic, Luipersbeck, Monaci, Resch, Salvagnin, Sinnl]
         */
        model_updated = separate_minimal_separators_integral(cuts_lhs, cuts_rhs);
    }
    else
    {
        /* run separation algorithm from "Partitioning a graph into balanced
         * connected classes - Formulations, separation and experiments", 2021,
         * by [Miyazawa, Moura, Ota, Wakabayashi]
         */
        model_updated = separate_minimal_separators_std(cuts_lhs, cuts_rhs);
    }

    if (model_updated)
    {
        // add cuts
        for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
        {
            ++minimal_separators_counter;

            if (kind_of_cut == ADD_USER_CUTS)
                addCut(cuts_lhs[idx] <= cuts_rhs[idx]);

            else if (kind_of_cut == ADD_LAZY_CNTRS)
                addLazy(cuts_lhs[idx] <= cuts_rhs[idx]);

            else // kind_of_cut == ADD_STD_CNTRS
                model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);
        }
    }

    return model_updated;
}

bool WCMCutGenerator::separate_minimal_separators_integral(vector<GRBLinExpr> &cuts_lhs,
                                                           vector<long> &cuts_rhs)
{
    /***
     * Solve the separation problem for minimal (a,b)-separator inequalities,
     * assuming the current point is integral
     */

    vector<long> vars_at_one = vector<long>();
    for (long u=0; u < num_vertices; ++u)
        if (y_val[u] >= MSI_ONE)
            vars_at_one.push_back(u);

    long num_vars_at_one = vars_at_one.size();
    if (num_vars_at_one < 2)
        return false;

    // 1. SUBGRAPH CONTAINING ONLY EDGES BETWEEN VERTICES AT ONE
    vector< vector<long> > aux_adj_list;

    // all vertices
    for (long i = 0; i < num_vertices; ++i)
        aux_adj_list.push_back(vector<long>());

    // only edges between vertices at one
    for (long i = 0; i < num_vars_at_one; ++i)
        for (long j = i+1; j < num_vars_at_one; ++j)
        {
            long u = vars_at_one.at(i);
            long v = vars_at_one.at(j);
            if (instance->graph->index_matrix[u][v] >= 0)
            {
                // i-th vertex at one (u) adjacent to j-th one (v)
                aux_adj_list[u].push_back(v);
                aux_adj_list[v].push_back(u);
            }
        }

    // 2. DFS IN THIS AUXILIARY GRAPH TAGGING CONNECTED COMPONENTS
    vector<long> components = vector<long>(num_vertices, -1);
    check_components(aux_adj_list, components);

    // 3. GET TWO VARS Y_s = Y_t = 1, WITH s AND t IN DIFFERENT COMPONENTS
    // NB! Trying to stick to the "rotating source" strategy to avoid favouring
    // separators between vertices of smaller index
    long s = -1;
    vector<long>::iterator it = vars_at_one.begin();
    while ( it != vars_at_one.end() && s < 0)
    {
        long v = *it;
        if(v >= this->msi_next_source)
            s = v;

        ++it;
    }

    // msi_next_source not at 1, nor any vertex with larger index?
    if (s < 0)
        s = vars_at_one.front();

    long t = -1;
    it = vars_at_one.begin();
    while ( it != vars_at_one.end() && t < 0)
    {
        long v = *it;
        if(components.at(v) != components.at(s))
            t = v;

        ++it;
    }

    if (t < 0)
    {
        // no two vertices at 1 in two different components... INTEGER FEASIBLE POINT!
        //cout << "### no two vertices at 1 in two different components... INTEGER FEASIBLE POINT" << endl;
        return false;
    }

    // 4. DETERMINE VERTICES IN V\COMPONENT[s] THAT ARE ADJACENT TO SOME VERTEX IN COMPONENT[s]
    vector<long> separator_vertices = vector<long>();
    vector<bool> separator_mask = vector<bool>(num_vertices, false);
    vector<bool> s_component_mask = vector<bool>(num_vertices, false);
    for (long u=0; u < num_vertices; ++u)
    {
        if(components.at(u) == components.at(s))
            s_component_mask.at(u) = true;
        else
        {
            bool u_is_a_neighbour = false;
            list<long>::iterator it = instance->graph->adj_list.at(u).begin();
            while (it != instance->graph->adj_list.at(u).end() && !u_is_a_neighbour)
            {
                long v = *it;
                if (components.at(v) == components.at(s))
                    u_is_a_neighbour = true;

                ++it;
            }

            if (u_is_a_neighbour)
            {
                separator_vertices.push_back(u);
                separator_mask.at(u) = true;
            }
        }
    }

    // 5. FOUND A SEPARATOR, BUT NOW LIFT IT TO A MINIMAL ONE
    #ifdef DEBUG_MSI_INTEGRAL
        cout << "### (" << s << "," << t << ")- separator"
             << endl;
        cout << "### before lifting: { ";

        for (vector<long>::iterator it = separator_vertices.begin();
                                    it != separator_vertices.end(); ++it)
            cout << *it << " ";

        cout << "}" << endl;
    #endif

    lift_to_minimal_separator(separator_vertices, separator_mask, s, t);

    #ifdef DEBUG_MSI_INTEGRAL
        cout << "### after lifting: { ";

        for (vector<long>::iterator it = separator_vertices.begin();
                                    it != separator_vertices.end(); ++it)
            cout << *it << " ";

        cout << "}" << endl;
    #endif

    // 6. DETERMINE INEQUALITY

    GRBLinExpr violated_constr = 0;

    violated_constr += y_vars[s];
    violated_constr += y_vars[t];

    vector<long>::iterator it_S = separator_vertices.begin();
    while (it_S != separator_vertices.end())
    {
        violated_constr += ( (-1) * y_vars[*it_S] );
        ++it_S;
    }

    cuts_lhs.push_back(violated_constr);
    cuts_rhs.push_back(1);

    #ifdef DEBUG_MSI_INTEGRAL
        double violating_lhs = 0;

        cout << "### ADDED MSI: ";
        cout << "y_" << s << " + y_" << t;

        violating_lhs += y_val[s];
        violating_lhs += y_val[t];

        it_S = separator_vertices.begin();
        while (it_S != separator_vertices.end())
        {
            cout << " - y_" << *it_S << "";
            violating_lhs -= y_val[*it_S];
            ++it_S;
        }

        cout << " <= 1 " << endl;
        cout << right;
        cout << setw(80) << "(lhs at current point "
             << violating_lhs << ")" << endl;
        cout << left;
    #endif

    this->msi_next_source++;
    return true;
}

bool WCMCutGenerator::separate_minimal_separators_std(vector<GRBLinExpr> &cuts_lhs,
                                                      vector<long> &cuts_rhs)
{
    /// Solve the separation problem for minimal (a,b)-separator inequalities

    // 1. CONSTRUCT AUXILIARY NETWORK D, WITH REDUCTIONS FROM INTEGRAL VARS

    // LEMON digraph representing the current solution
    SmartDigraph D;
    vector<SmartDigraph::Node> D_vertices;
    map<pair<long,long>, SmartDigraph::Arc> D_arcs;
    SmartDigraph::ArcMap<double> D_capacity(D);
    long D_size = 0;

    vector<long> vars_at_one = vector<long>();
    vector<long> fractional_vars_D_idx = vector<long>();
    vector<double> fractional_vars_val = vector<double>();

    // maps u->u_1 ; u_2 = D_idx_of_vertex[u]+1 for fractional y_val[u]
    vector<long> D_idx_of_vertex = vector<long>(num_vertices, -1);

    // 1.1 ADD VERTICES CORRESPONDING TO VARS IN [0,1) IN THIS RELAXATION
    for (long u = 0; u < num_vertices; ++u)
    {
        const double value = y_val[u];
        if (value <= MSI_ZERO)  // == 0
        {
            // add only one vertex u_1 = u_2 in D
            D_vertices.push_back(D.addNode());
            D_idx_of_vertex[u] = D_size;
            ++D_size;
        }
        else if(value > MSI_ZERO && value < MSI_ONE)   // > 0 && < 1
        {
            // add two vertices u_1, u_2 in D
            D_vertices.push_back(D.addNode());
            D_vertices.push_back(D.addNode());
            D_idx_of_vertex[u] = D_size;

            // redundant, but helps adding arcs more efficiently below
            fractional_vars_D_idx.push_back(D_size);
            fractional_vars_val.push_back(value);

            D_size += 2;
        }
        else // value == 1
            vars_at_one.push_back(u);
    }

    // 1.2 ADD VERTICES CORRESPONDING TO VARS AT 1

    // inspect subgraph induced by vertices at one (contracted if adjacent)
    long num_vars_at_one = vars_at_one.size();

    vector< vector<long> > aux_adj_list;

    for (long i = 0; i < num_vars_at_one; ++i)
        aux_adj_list.push_back(vector<long>());

    for (long i = 0; i < num_vars_at_one; ++i)
        for (long j = i+1; j < num_vars_at_one; ++j)
        {
            long u = vars_at_one.at(i);
            long v = vars_at_one.at(j);
            if (instance->graph->index_matrix[u][v] >= 0)
            {
                // i-th vertex at one (u) adjacent to j-th one (v)
                aux_adj_list[i].push_back(j);
                aux_adj_list[j].push_back(i);
            }
        }

    // dfs in this auxiliary graph tagging connected components
    vector<long> components = vector<long>(num_vars_at_one, -1);
    long num_components = check_components(aux_adj_list, components);

    // add one vertex in D for each component in the auxiliary graph
    for (long i = 0; i < num_components; ++i)
        D_vertices.push_back(D.addNode());

    for (long i = 0; i < num_vars_at_one; ++i)
    {
        long u = vars_at_one.at(i);
        long cluster = D_size + components.at(i); // D_size not updated yet
        D_idx_of_vertex[u] = cluster;
    }

    D_size += num_components;

    // 1.3 FIRST SET OF ARCS: (U1,U2) FOR U IN V(G) S.T. y_val[u] \in (0,1)

    const unsigned num_frac_vars = fractional_vars_D_idx.size();

    for (unsigned idx = 0; idx < num_frac_vars; ++idx)
    {
        long v1 = fractional_vars_D_idx.at(idx);
        long v2 = v1 + 1;
        pair<long,long> v1v2 = make_pair(v1,v2);

        D_arcs[v1v2] = D.addArc(D_vertices[v1], D_vertices[v2]);

        D_capacity[ D_arcs[v1v2] ] = fractional_vars_val.at(idx);
    }

    // 1.4 REMAINING ARCS: UNLIMITED CAPACITY (|V|+1 SUFFICES HERE!)

    const long num_edges = instance->graph->num_edges;
    const long UNLIMITED_CAPACITY = num_vertices + 1;

    for (long idx = 0; idx < num_edges; ++idx)
    {
        long u = instance->graph->s.at(idx);
        long v = instance->graph->t.at(idx);

        // the actual arcs (one or two, for each original edge) depend on
        // the reductions due to integral valued vars (7 cases... boring!)
        // NB!
        // "... >= MSI_ONE" means ... == 1; "... <= MSI_ZERO" means == 0
        // "... >= MSI_ONE" means ... == 1; "... <= MSI_ZERO" means == 0
        // "... >= MSI_ONE" means ... == 1; "... <= MSI_ZERO" means == 0

        if (y_val[u] >= MSI_ONE && y_val[v] <= MSI_ZERO)
        {
            // CASE 1
            long uu = D_idx_of_vertex[u];
            long vv = D_idx_of_vertex[v];

            pair<long,long> uuvv = make_pair(uu,vv);
            D_arcs[uuvv] = D.addArc(D_vertices[uu], D_vertices[vv]);
            D_capacity[ D_arcs[uuvv] ] = UNLIMITED_CAPACITY;

            //pair<long,long> vvuu = make_pair(vv,uu);
            //D_arcs[vvuu] = D.addArc(D_vertices[vv], D_vertices[uu]);
            //D_capacity[ D_arcs[vvuu] ] = UNLIMITED_CAPACITY;
        }
        else if (y_val[u] <= MSI_ZERO && y_val[v] >= MSI_ONE)
        {
            // CASE 2
            long uu = D_idx_of_vertex[u];
            long vv = D_idx_of_vertex[v];

            pair<long,long> vvuu = make_pair(vv,uu);
            D_arcs[vvuu] = D.addArc(D_vertices[vv], D_vertices[uu]);
            D_capacity[ D_arcs[vvuu] ] = UNLIMITED_CAPACITY;

            //pair<long,long> uuvv = make_pair(uu,vv);
            //D_arcs[uuvv] = D.addArc(D_vertices[uu], D_vertices[vv]);
            //D_capacity[ D_arcs[uuvv] ] = UNLIMITED_CAPACITY;
        }
        else if (y_val[u] <= MSI_ZERO && y_val[v] > MSI_ZERO
                                      && y_val[v] < MSI_ONE)
        {
            // CASE 3
            long uu = D_idx_of_vertex[u];
            //long v1 = D_idx_of_vertex[v];
            long v2 = D_idx_of_vertex[v] + 1;

            pair<long,long> v2uu = make_pair(v2,uu);
            D_arcs[v2uu] = D.addArc(D_vertices[v2], D_vertices[uu]);
            D_capacity[ D_arcs[v2uu] ] = UNLIMITED_CAPACITY;

            //pair<long,long> uuv1 = make_pair(uu,v1);
            //D_arcs[uuv1] = D.addArc(D_vertices[uu], D_vertices[v1]);
            //D_capacity[ D_arcs[uuv1] ] = UNLIMITED_CAPACITY;
        }
        else if (y_val[u] > MSI_ZERO && y_val[v] <= MSI_ZERO &&
                 y_val[u] < MSI_ONE)
        {
            // CASE 4
            //long u1 = D_idx_of_vertex[u];
            long u2 = D_idx_of_vertex[u] + 1;
            long vv = D_idx_of_vertex[v];

            pair<long,long> u2vv = make_pair(u2,vv);
            D_arcs[u2vv] = D.addArc(D_vertices[u2], D_vertices[vv]);
            D_capacity[ D_arcs[u2vv] ] = UNLIMITED_CAPACITY;

            //pair<long,long> vvu1 = make_pair(vv,u1);
            //D_arcs[vvu1] = D.addArc(D_vertices[vv], D_vertices[u1]);
            //D_capacity[ D_arcs[vvu1] ] = UNLIMITED_CAPACITY;
        }
        else if (y_val[u] > MSI_ZERO && y_val[v] >= MSI_ONE &&
                 y_val[u] < MSI_ONE)
        {
            // CASE 5
            long u1 = D_idx_of_vertex[u];
            long u2 = D_idx_of_vertex[u] + 1;
            long vv = D_idx_of_vertex[v];
            pair<long,long> u2vv = make_pair(u2,vv);
            pair<long,long> vvu1 = make_pair(vv,u1);

            D_arcs[u2vv] = D.addArc(D_vertices[u2], D_vertices[vv]);
            D_arcs[vvu1] = D.addArc(D_vertices[vv], D_vertices[u1]);

            D_capacity[ D_arcs[u2vv] ] = UNLIMITED_CAPACITY;
            D_capacity[ D_arcs[vvu1] ] = UNLIMITED_CAPACITY;
        }
        else if (y_val[u] >= MSI_ONE && y_val[v] > MSI_ZERO
                                     && y_val[v] < MSI_ONE)
        {
            // CASE 6
            long uu = D_idx_of_vertex[u];
            long v1 = D_idx_of_vertex[v];
            long v2 = D_idx_of_vertex[v] + 1;
            pair<long,long> v2uu = make_pair(v2,uu);
            pair<long,long> uuv1 = make_pair(uu,v1);

            D_arcs[v2uu] = D.addArc(D_vertices[v2], D_vertices[uu]);
            D_arcs[uuv1] = D.addArc(D_vertices[uu], D_vertices[v1]);

            D_capacity[ D_arcs[v2uu] ] = UNLIMITED_CAPACITY;
            D_capacity[ D_arcs[uuv1] ] = UNLIMITED_CAPACITY;
        }
        else if (y_val[u] > MSI_ZERO && y_val[v] > MSI_ZERO &&
                 y_val[u] < MSI_ONE  && y_val[v] < MSI_ONE  )
        {
            // CASE 7
            long u1 = D_idx_of_vertex[u];
            long u2 = D_idx_of_vertex[u] + 1;
            long v1 = D_idx_of_vertex[v];
            long v2 = D_idx_of_vertex[v] + 1;
            pair<long,long> u2v1 = make_pair(u2,v1);
            pair<long,long> v2u1 = make_pair(v2,u1);

            D_arcs[u2v1] = D.addArc(D_vertices[u2], D_vertices[v1]);
            D_arcs[v2u1] = D.addArc(D_vertices[v2], D_vertices[u1]);

            D_capacity[ D_arcs[u2v1] ] = UNLIMITED_CAPACITY;
            D_capacity[ D_arcs[v2u1] ] = UNLIMITED_CAPACITY;
        }
    }

    /* 2. TRY EACH PAIR OF NON-ADJACENT VERTICES IN THE ORIGINAL GRAPH WHOSE
     * COMBINED Y VALUES IN THIS RELAXATION SOLUTION EXCEED 1 (THE CORRESPONDING 
     * MSI CANNOT BE VIOLATED OTHERWISE)
     */

    bool done = false;
    long num_trials = 0;
    while (num_trials < num_vertices && !done)
    {
        // "cycling" through the initial source vertex tried (only relevant with MSI_STRATEGY_FIRST_CUT_BELOW_ROOT)
        long s = this->msi_next_source;
        long t = s+1;
        while (t < num_vertices && !done)
        {
            // wanted: a (s_2, t_1) separating cut
            long s_in_D = (y_val[s] > MSI_ZERO && y_val[s] < MSI_ONE) ? D_idx_of_vertex[s]+1
                                                                      : D_idx_of_vertex[s];

            long t_in_D = D_idx_of_vertex[t];

            if ( instance->graph->index_matrix[s][t] < 0 &&  // non-adjacent
                 y_val[s] + y_val[t] > 1+MSI_EPSILON &&      // might cut y*
                 s_in_D != t_in_D )                          // not contracted
            {
                // 3. MAX FLOW COMPUTATION

                /***
                 * Using the first phase of Goldberg & Tarjan preflow
                 * push-relabel algorithm (with "highest label" and "bound
                 * decrease" heuristics). The worst case time complexity of
                 * the algorithm is in O(n^2 * m^0.5), n and m wrt D
                 */

                Preflow<SmartDigraph, SmartDigraph::ArcMap<double> >
                    s_t_preflow(D, D_capacity, D_vertices[s_in_D],
                                               D_vertices[t_in_D]);

                s_t_preflow.runMinCut();
                
                double mincut = s_t_preflow.flowValue();

                // 4. IF THE MIN CUT IS LESS THAN WHAT THE MSI PRESCRIBES
                // (UP TO A VIOLATION TOLERANCE), WE FOUND A CUT

                if (mincut < y_val[s] + y_val[t] - 1 - MSI_EPSILON)
                {
                    // 5. DETERMINE VERTICES IN ORIGINAL GRAPH CORRESPONDING
                    // TO ARCS IN THE MIN CUT

                    vector<long> S = vector<long>();
                    vector<bool> S_mask = vector<bool>(num_vertices, false);

                    for (long u = 0; u < num_vertices; ++u)
                    {
                        if (u != s && u != t)
                        {
                            long u_in_D = D_idx_of_vertex[u];

                            // query if u1 is on the source side of the min cut
                            if ( s_t_preflow.minCut(D_vertices[u_in_D]) )
                            {
                                bool u_at_zero = y_val[u] <= MSI_ZERO;

                                bool u_frac = y_val[u] > MSI_ZERO && y_val[u] < MSI_ONE;

                                long u2_in_D = u_frac ? D_idx_of_vertex[u]+1
                                                      : D_idx_of_vertex[u];

                                bool u2_separated =
                                   !s_t_preflow.minCut(D_vertices[u2_in_D]);

                                if (u_at_zero || (u_frac && u2_separated) )
                                {
                                    S.push_back(u);
                                    S_mask.at(u) = true;
                                }
                            }
                        }
                    }

                    // 6. LIFT CUT BY REDUCING S TO A MINIMAL SEPARATOR
                    #ifdef DEBUG_MSI
                        cout << "### (" << s << "," << t << ")- separator"
                             << endl;
                        cout << "### before lifting: { ";

                    for (vector<long>::iterator it = S.begin();
                                                it != S.end(); ++it)
                        cout << *it << " ";

                    cout << "}" << endl;
                    #endif

                    lift_to_minimal_separator(S, S_mask, s, t);

                    #ifdef DEBUG_MSI
                        cout << "### after lifting: { ";

                    for (vector<long>::iterator it = S.begin();
                                                it != S.end(); ++it)
                        cout << *it << " ";

                    cout << "}" << endl;
                    #endif

                    // 7. DETERMINE INEQUALITY

                    GRBLinExpr violated_constr = 0;

                    violated_constr += y_vars[s];
                    violated_constr += y_vars[t];

                    for (vector<long>::iterator it = S.begin();
                                                it != S.end(); ++it)
                    {
                        violated_constr += ( (-1) * y_vars[*it] );
                    }

                    cuts_lhs.push_back(violated_constr);
                    cuts_rhs.push_back(1);

                    #ifdef DEBUG_MSI
                        double violating_lhs = 0;

                        cout << "### ADDED MSI: ";
                        cout << "y_" << s << " + y_" << t;

                        violating_lhs += y_val[s];
                        violating_lhs += y_val[t];

                        for (vector<long>::iterator it = S.begin();
                                                    it != S.end(); ++it)
                        {
                            cout << " - y_" << *it << "";
                            violating_lhs -= y_val[*it];
                        }

                        cout << " <= 1 " << endl;
                        cout << right;
                        cout << setw(80) << "(lhs at current point "
                             << violating_lhs << ")" << endl;
                        cout << left;
                    #endif

                    if (MSI_STRATEGY_FIRST_CUT_BELOW_ROOT && !at_root_relaxation)
                        done = true;
                }
            }

            ++t;
        }

        this->msi_next_source = (this->msi_next_source+1) % num_vertices;
        ++num_trials;
    }

    return (cuts_lhs.size() > 0);
}

void inline WCMCutGenerator::lift_to_minimal_separator(vector<long> &S,
                                                       vector<bool> &S_mask,
                                                       long s,
                                                       long t)
{
    /// Remove vertices from S until it is a minimal (s,t)-separator

    bool updated = true;
    while (updated)
    {
        updated = false;

        const long len = S.size();

        long count = 0;
        vector<bool> seen = vector<bool>(num_vertices, false);

        // only try dfs from t if dfs from s found all vertices of S 
        dfs_avoiding_set(S, S_mask, s, seen, count);
        if (count == len)
        {
            count = 0;
            seen.assign(num_vertices, false);
            dfs_avoiding_set(S, S_mask, t, seen, count);
        }

        /***
         * S is a minimal (s,t)-separator if and only if every element in S has
         * a neighbour both in the connected component of G-S containing s, and  
         * in the one containing t. So we may remove from S a vertex not found
         * in either dfs' above.
         */
        if (count < len)
        {
            long i = 0;
            while (i < len && !updated)
            {
                long vertex_at_i = S.at(i);

                if ( !seen.at(vertex_at_i) )
                {
                    updated = true;
                    S.erase(S.begin() + i);
                    S_mask.at(vertex_at_i) = false;
                }

                ++i;
            }
        }

    }  // repeat search if S was updated 
}

void inline WCMCutGenerator::dfs_avoiding_set(vector<long> &S,
                                              vector<bool> &S_mask,
                                              long source,
                                              vector<bool> &seen,
                                              long &count)
{
    // auxiliary dfs tagging seen vertices, but not exploring vertices in S

    // NB! Expecting 'seen' with num_vertices 'false' entries (not S.size()!)

    seen.at(source) = true;

    for (list<long>::iterator it = instance->graph->adj_list.at(source).begin();
        it != instance->graph->adj_list.at(source).end(); ++it)
    {
        long v = *it;

        if ( !seen.at(v) )
        {
            if ( S_mask.at(v) )
            {
                // v in S
                seen.at(v) = true;
                ++count;
            }
            else
                dfs_avoiding_set(S, S_mask, v, seen, count);
        }
    }
}

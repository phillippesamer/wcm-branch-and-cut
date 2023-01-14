#include "cks_cutgenerator.h"

/// algorithm setup switches

bool SEPARATE_MSI = true;               // MSI = minimal separator inequalities
bool SEPARATE_INDEGREE = true;

bool CUTS_AT_ROOT_ONLY = false;

// strategy for running separation algorithms for colour-specific inequalities
bool SEARCH_ALL_COLOURS_FOR_INDEGREE = true;
bool SEARCH_ALL_COLOURS_FOR_MSI = true;

// at most 14 without changing everything to long double (which gurobi ignores)
#define SEPARATION_PRECISION_IN_IP 14

///////////////////////////////////////////////////////////////////////////////

/// specialized depth-first search to identify/count connected components

void dfs_to_tag_components(long u,
                           long count, 
                           vector<long> &components, 
                           vector< vector<long> > &adj_list)
{
    // auxiliary dfs to check connected components in adj_list
    components[u] = count;

    for (unsigned i=0; i<adj_list[u].size(); ++i)
    {
        long v = adj_list[u].at(i);
        if (components[v] < 0)
            dfs_to_tag_components(v, count, components, adj_list);
    }
}

long check_components(vector< vector<long> > &adj_list,
                      vector<long> &components)
{
    long count = 0;
    for (long u = 0; u < adj_list.size(); ++u)
    {
        if (components[u] < 0)
            dfs_to_tag_components(u, count, components, adj_list);

        ++count;
    }

    return count;
}

///////////////////////////////////////////////////////////////////////////////

CKSCutGenerator::CKSCutGenerator(GRBModel *model, GRBVar **x_vars, IO *instance)
{
    this->model = model;
    this->x_vars = x_vars;
    this->instance = instance;

    this->num_vertices = instance->graph->num_vertices;
    this->num_subgraphs = instance->num_subgraphs;

    this->indegree_counter = 0;
    this->minimal_separators_counter = 0;
}

void CKSCutGenerator::callback()
{
    /***
     * The actual callback method within the solver. Currently, only used for 
     * dynamically adding violated indegree inequalities (if any) and/or
     * violated minimal separator inequalities (MSI).
     * NB! SEARCHING FOR MSI ONLY WHEN NO VIOLATED INDEGREE EXISTS.
     */

    try
    {
        // callback from the search at a given MIP node: including USER CUTS
        if (where == GRB_CB_MIPNODE)
        {
            // node relaxation solution must be available at the current node
            if (this->getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL)
                return;

            // generate cuts only at root node?
            if (CUTS_AT_ROOT_ONLY && getDoubleInfo(GRB_CB_MIPNODE_NODCNT) > 0)
                return;

            // retrieve relaxation solution
            x_val = new double*[num_vertices];
            for (long u = 0; u < num_vertices; ++u)
                x_val[u] = this->getNodeRel(x_vars[u], num_subgraphs);

            clean_x_val_beyond_precision(SEPARATION_PRECISION_IN_IP);

            bool model_updated = false;

            if (SEPARATE_INDEGREE)
                model_updated = run_indegree_separation(ADD_USER_CUTS);

            if (!model_updated && SEPARATE_MSI)
                run_minimal_separators_separation(ADD_USER_CUTS);

            for (long u=0; u < num_vertices; u++)
                delete[] x_val[u];
            delete[] x_val;
        }

        // callback from a new MIP incumbent: including LAZY CONSTRAINTS
        else if (where == GRB_CB_MIPSOL)
        {
            // retrieve solution
            x_val = new double*[num_vertices];
            for (long u = 0; u < num_vertices; ++u)
                x_val[u] = this->getSolution(x_vars[u], num_subgraphs);

            clean_x_val_beyond_precision(SEPARATION_PRECISION_IN_IP);

            bool model_updated = false;

            if (SEPARATE_INDEGREE)
                model_updated = run_indegree_separation(ADD_LAZY_CNTRS);

            if (!model_updated && SEPARATE_MSI)
                run_minimal_separators_separation(ADD_LAZY_CNTRS);

            for (long u=0; u < num_vertices; u++)
                delete[] x_val[u];
            delete[] x_val;
        }
    }
    catch (GRBException e)
    {
        cout << "Error " << e.getErrorCode()
             << " during CKSCutGenerator::callback(): ";
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Unexpected error during CKSCutGenerator::callback()" << endl;
    }
}

bool CKSCutGenerator::separate_lpr()
{
    /// Interface to be used when solving the LP relaxation only.

    try
    {
        // retrieve relaxation solution
        x_val = new double*[num_vertices];
        for (long u = 0; u < num_vertices; ++u)
        {
            x_val[u] = new double[num_subgraphs];
            for (long c = 0; c < num_subgraphs; ++c)
                x_val[u][c] = x_vars[u][c].get(GRB_DoubleAttr_X);
        }

        clean_x_val_beyond_precision(SEPARATION_PRECISION_IN_IP);

        bool model_updated = false;

        if (SEPARATE_INDEGREE)
            model_updated = run_indegree_separation(ADD_STD_CNTRS);

        if (!model_updated && SEPARATE_MSI)
            model_updated = run_minimal_separators_separation(ADD_STD_CNTRS);

        for (long u=0; u < num_vertices; u++)
            delete[] x_val[u];
        delete[] x_val;

        return model_updated;
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

void CKSCutGenerator::clean_x_val_beyond_precision(int precision)
{
    /// prevent floating point errors by ignoring digits beyond given precision
    for (long u = 0; u < num_vertices; ++u)
        for (long c = 0; c < num_subgraphs; ++c)
        {
            double tmp = x_val[u][c] * std::pow(10, precision);
            tmp = std::round(tmp);
            x_val[u][c] = tmp * std::pow(10, -precision);
        }
}

bool CKSCutGenerator::run_indegree_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    /* run separation algorithm from "On imposing connectivity constraints in
     * integer programs", 2017, by [Wang, Buchanan, Butenko]
     */
    model_updated = separate_indegree(cuts_lhs,cuts_rhs);

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

bool CKSCutGenerator::separate_indegree(vector<GRBLinExpr> &cuts_lhs,
                                        vector<long> &cuts_rhs)
{
    /// Solve the separation problem for indegree inequalities, for each colour

    const long num_edges = instance->graph->num_edges;

    long colour = 0;
    bool done = false;
    while (colour < num_subgraphs && !done)
    {
        vector<long> indegree = vector<long>(num_vertices, 0);

        // 1. COMPUTE INDEGREE ORIENTING EDGES ACCORDING TO RELAXATION SOLUTION
        for (long idx = 0; idx < num_edges; ++idx)
        {
            long u = instance->graph->s.at(idx);
            long v = instance->graph->t.at(idx);
            if (x_val[u][colour] > x_val[v][colour])
                indegree.at(v) += 1;
            else
                indegree.at(u) += 1;
        }

        // 2. EVALUATE LHS (1-d[u])*x[u]
        double lhs_sum = 0.;
        for (long u = 0; u < num_vertices; ++u)
            lhs_sum += ( (1 - indegree.at(u)) * x_val[u][colour] );

        // 3. FOUND MOST VIOLATED INDEGREE INEQUALITY (IF ANY) IF LHS > 1
        if (lhs_sum > 1)
        {
            // store inequality (caller method adds it to the model)
            GRBLinExpr violated_constr = 0;

            for (long u = 0; u < num_vertices; ++u)
                violated_constr += ( (1 - indegree.at(u)) * x_vars[u][colour] );

            cuts_lhs.push_back(violated_constr);
            cuts_rhs.push_back(1);

            if (!SEARCH_ALL_COLOURS_FOR_INDEGREE)
                done = true;
        }

        ++colour;
    }

    return (cuts_lhs.size() > 0);
}

bool CKSCutGenerator::run_minimal_separators_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    /* run separation algorithm from "Partitioning a graph into balanced
     * connected classes - Formulations, separation and experiments", 2021,
     * by [Miyazawa, Moura, Ota, Wakabayashi]
     */
    model_updated = separate_minimal_separators(cuts_lhs,cuts_rhs);

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

bool CKSCutGenerator::separate_minimal_separators(vector<GRBLinExpr> &cuts_lhs,
                                                  vector<long> &cuts_rhs)
{
    /// Solve the separation problem for minimal (a,b)-separator inequalities

    long colour = 0;
    bool done = false;
    while (colour < num_subgraphs && !done)
    {
        // 1. CONSTRUCT AUXILIARY NETWORK D, WITH REDUCTIONS FROM INTEGRAL VARS

        // LEMON digraph representing the current solution
        SmartDigraph lemon_g;
        SmartDigraph::ArcMap<double> lemon_cap(lemon_g);
        vector<SmartDigraph::Node> lemon_vertices;
        map<pair<long,long>, SmartDigraph::Arc> lemon_arcs;
        long size_D = 0;

        vector<long> vars_at_one = vector<long>();

        // maps u->u_1 ; u_2 = D_idx_of_vertex[u]+1 for fractional x_val[u]
        vector<long> D_idx_of_vertex = vector<long>(num_vertices, -1);

        for (long u = 0; u < num_vertices; ++u)
        {
            if (x_val[u][colour] == 0)
            {
                // add only one vertex u_1 = u_2 in D
                lemon_vertices.push_back(lemon_g.addNode());
                D_idx_of_vertex[u] = size_D;
                ++size_D;
            }
            else if(x_val[u][colour] > 0 && x_val[u][colour] < 1)
            {
                // add two vertices u_1, u_2 in D
                lemon_vertices.push_back(lemon_g.addNode());
                lemon_vertices.push_back(lemon_g.addNode());
                D_idx_of_vertex[u] = size_D;
                size_D += 2;
            }
            else // x_val[u][colour] == 1
                vars_at_one.push_back(u);
        }

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
        vector<long> components = vector<long>(num_var_at_one, -1);
        long num_components = check_components(aux_adj_list, components);

        // add one vertex in D for each component in the auxiliary graph
        for (long i = 0; i < num_components; ++i)
            lemon_vertices.push_back(lemon_g.addNode());

        for (long i = 0; i < num_vars_at_one; ++i)
        {
            long u = vars_at_one.at(i);
            long cluster = size_D + components.at(i); // size_D not updated yet
            D_idx_of_vertex[u] = cluster;
        }

        size_D += num_components;

        // TO DO: check implementation of check_components
        // !!! CUIDADO COM ÍNDICES i,j ACIMA E ÍNDICES REAIS (ARMAZENADOS EM vars_at_one)






        // create arcs



        // lemon_arcs[make_pair(v1,v2)] = lemon_g.addArc(lemon_vertices[v1], lemon_vertices[v2]);
        // lemon_cap[ lemon_arcs[make_pair(v1,v2)] ] = x_val[v][colour];





        // TO DO: replace names of lemon structures




        /* 2. TRY EACH PAIR OF NON-ADJACENT VERTICES (IN THE ORIGINAL GRAPH)
         * WHOSE COMBINED VALUES IN THIS RELAXATION SOLUTION EXCEED 1
         * AS LONG AS A VIOLATED MSI IS NOT FOUND
         */

        bool done_with_this_colour = false;
        long s = 0;
        while (s < num_vertices && !done_with_this_colour)
        {
            long t = s+1;
            while (t < num_vertices && !done_with_this_colour)
            {
                if ( instance->graph->index_matrix[s][t] < 0 &&
                     x_val[s][colour] + x_val[t][colour] > 1 )
                {
                    // 3. MAX FLOW COMPUTATION



                    // 4. IF THE MAX FLOW (MIN CUT) IS LESS THAN WHAT THE
                    // INEQUALITY PRESCRIBES, WE FOUND A CUT

                    if (maxflow_value < x_val[s][colour] + x_val[t][colour] - 1)
                    {
                        // 5. DETERMINE VERTICES IN ORIGINAL GRAPH CORRESPONDING
                        // TO ARCS IN THE MIN CUT




                        // TO DO: 6. LIFT VIOLATED INEQUALITY BY REDUCING S






                        // determine inequality







                        done_with_this_colour = true;

                        //if (!SEARCH_ALL_COLOURS_FOR_MSI)
                        //    done = true;
                    }





                }







                ++t;
            }

            ++s;
        }







        ++colour;
    }

    return (cuts_lhs.size() > 0);
}

#include "cks_cutgenerator.h"

/// algorithm setup switches

bool SEPARATE_MSI = true;               // MSI = minimal separator inequalities
bool SEPARATE_INDEGREE = false;

bool MSI_ONLY_IF_NO_INDEGREE = false;   // only used with SEPARATE_MSI == true

// strategy for running separation algorithms for colour-specific inequalities
bool SEARCH_ALL_COLOURS_FOR_INDEGREE = true;
bool SEARCH_ALL_COLOURS_FOR_MSI = false;

// at most 14 without changing everything to long double (which gurobi ignores)
bool SET_MAX_PRECISION_IN_SEPARATION = true;
int  SEPARATION_PRECISION = 14;
double MSI_VIOLATION_TOL = 1e-10;

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

    for (long i = 0; i < degree; ++i)
    {
        long v = adj_list[u].at(i);
        if (components[v] < 0)
            dfs_to_tag_components(v, count, components, adj_list);
    }
}

long inline check_components(vector< vector<long> > &adj_list,
                             vector<long> &components)
{
    /// NB! Expect components initialized as vector<long>(adj_list.size(), -1)

    long count = 0;
    const long n = adj_list.size();

    for (long u = 0; u < n; ++u)
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
    this->at_root_relaxation = true;

    this->indegree_counter = 0;
    this->minimal_separators_counter = 0;

    this->msi_current_colour = 0;
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

            // flag when done with the root node relaxation
            if (this->at_root_relaxation)   // initially true
                if (getDoubleInfo(GRB_CB_MIPNODE_NODCNT) > 0)
                    this->at_root_relaxation = false;

            // retrieve relaxation solution
            x_val = new double*[num_vertices];
            for (long u = 0; u < num_vertices; ++u)
                x_val[u] = this->getNodeRel(x_vars[u], num_subgraphs);

            bool separated = false;

            if (SEPARATE_INDEGREE)
                separated = run_indegree_separation(ADD_USER_CUTS);

            if (SEPARATE_MSI)
            {
                if (!MSI_ONLY_IF_NO_INDEGREE || !separated)
                {
                    if (SET_MAX_PRECISION_IN_SEPARATION)
                        clean_x_val_beyond_precision(SEPARATION_PRECISION);

                    run_minimal_separators_separation(ADD_USER_CUTS);
                }
            }

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

            bool separated = false;

            //if (SEPARATE_INDEGREE)
            //    separated = run_indegree_separation(ADD_LAZY_CNTRS);

            if (SEPARATE_MSI)
            {
                if (!MSI_ONLY_IF_NO_INDEGREE || !separated)
                {
                    if (SET_MAX_PRECISION_IN_SEPARATION)
                        clean_x_val_beyond_precision(SEPARATION_PRECISION);

                    run_minimal_separators_separation(ADD_LAZY_CNTRS);
                }
            }

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

        bool indegree_cut = false;
        bool msi_cut = false;

        if (SEPARATE_INDEGREE)
            indegree_cut = run_indegree_separation(ADD_STD_CNTRS);

        if (SEPARATE_MSI)
        {
            if (!MSI_ONLY_IF_NO_INDEGREE || !indegree_cut)
            {
                if (SET_MAX_PRECISION_IN_SEPARATION)
                    clean_x_val_beyond_precision(SEPARATION_PRECISION);

                msi_cut = run_minimal_separators_separation(ADD_STD_CNTRS);
            }
        }

        for (long u=0; u < num_vertices; u++)
            delete[] x_val[u];
        delete[] x_val;

        return (indegree_cut || msi_cut);
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
            if (x_val[u][colour] > x_val[v][colour] + EPSILON_TOL)
                indegree.at(v) += 1;
            else
                indegree.at(u) += 1;
        }

        // 2. EVALUATE LHS (1-d[u])*x[u]
        double lhs_sum = 0.;
        for (long u = 0; u < num_vertices; ++u)
            lhs_sum += ( (1 - indegree.at(u)) * x_val[u][colour] );

        // 3. FOUND MOST VIOLATED INDEGREE INEQUALITY (IF ANY) IF LHS > 1
        if (lhs_sum > 1 + EPSILON_TOL)
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

    long colours_tried = 0;
    bool done = false;
    while (colours_tried < num_subgraphs && !done)
    {
        long colour = this->msi_current_colour;

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

        // maps u->u_1 ; u_2 = D_idx_of_vertex[u]+1 for fractional x_val[u]
        vector<long> D_idx_of_vertex = vector<long>(num_vertices, -1);

        // 1.1 ADD VERTICES CORRESPONDING TO VARS IN [0,1) IN THIS RELAXATION
        for (long u = 0; u < num_vertices; ++u)
        {
            const double value = x_val[u][colour];
            if (value == 0)
            {
                // add only one vertex u_1 = u_2 in D
                D_vertices.push_back(D.addNode());
                D_idx_of_vertex[u] = D_size;
                ++D_size;
            }
            else if(value > 0 && value < 1)
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

        // 1.3 FIRST SET OF ARCS: (U1,U2) FOR U IN V(G) S.T. x_val[u] \in (0,1)

        const unsigned num_frac_vars = fractional_vars_D_idx.size();

        for (unsigned idx = 0; idx < num_frac_vars; ++idx)
        {
            long v1 = fractional_vars_D_idx.at(idx);
            long v2 = v1 + 1;
            pair<long,long> v1v2 = make_pair(v1,v2);

            D_arcs[v1v2] = D.addArc(D_vertices[v1], D_vertices[v2]);

            D_capacity[ D_arcs[v1v2] ] = fractional_vars_val.at(idx);
        }

        // 1.4 REMAINING ARCS, WITH UNLIMITED CAPACITY (|V|+1 SUFFICES HERE!)

        const long num_edges = instance->graph->num_edges;
        const long UNLIMITED_CAPACITY = num_vertices + 1;

        for (long idx = 0; idx < num_edges; ++idx)
        {
            long u = instance->graph->s.at(idx);
            long v = instance->graph->t.at(idx);

            // the actual arcs (one or two, for each original edge) depend on
            // the reductions due to integral valued vars (7 cases... boring!)
            if (x_val[u][colour] == 1 && x_val[v][colour] == 0)
            {
                // CASE 1
                long uu = D_idx_of_vertex[u];
                long vv = D_idx_of_vertex[v];
                pair<long,long> uuvv = make_pair(uu,vv);

                D_arcs[uuvv] = D.addArc(D_vertices[uu], D_vertices[vv]);

                D_capacity[ D_arcs[uuvv] ] = UNLIMITED_CAPACITY;
            }
            else if (x_val[u][colour] == 0 && x_val[v][colour] == 1)
            {
                // CASE 2
                long uu = D_idx_of_vertex[u];
                long vv = D_idx_of_vertex[v];
                pair<long,long> vvuu = make_pair(vv,uu);

                D_arcs[vvuu] = D.addArc(D_vertices[vv], D_vertices[uu]);

                D_capacity[ D_arcs[vvuu] ] = UNLIMITED_CAPACITY;
            }
            else if (x_val[u][colour] == 0 && x_val[v][colour] > 0
                                           && x_val[v][colour] < 1)
            {
                // CASE 3
                long uu = D_idx_of_vertex[u];
                long v2 = D_idx_of_vertex[v] + 1;
                pair<long,long> v2uu = make_pair(v2,uu);

                D_arcs[v2uu] = D.addArc(D_vertices[v2], D_vertices[uu]);

                D_capacity[ D_arcs[v2uu] ] = UNLIMITED_CAPACITY;
            }
            else if (x_val[u][colour] > 0 && x_val[v][colour] == 0 &&
                     x_val[u][colour] < 1)
            {
                // CASE 4
                long u2 = D_idx_of_vertex[u] + 1;
                long vv = D_idx_of_vertex[v];
                pair<long,long> u2vv = make_pair(u2,vv);

                D_arcs[u2vv] = D.addArc(D_vertices[u2], D_vertices[vv]);

                D_capacity[ D_arcs[u2vv] ] = UNLIMITED_CAPACITY;
            }
            else if (x_val[u][colour] > 0 && x_val[v][colour] == 1 &&
                     x_val[u][colour] < 1)
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
            else if (x_val[u][colour] == 1 && x_val[v][colour] > 0
                                           && x_val[v][colour] < 1)
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
            else if (x_val[u][colour] > 0 && x_val[v][colour] > 0 &&
                     x_val[u][colour] < 1 && x_val[v][colour] < 1  )
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

        /* 2. TRY EACH PAIR OF NON-ADJACENT VERTICES (IN THE ORIGINAL GRAPH)
         * WHOSE COMBINED VALUES IN THIS RELAXATION SOLUTION EXCEED 1
         * (OTHERWISE THE CORRESPONDNG MSI CANNOT BE VIOLATED)
         */

        // NOTE: not looking for more than 1 violated MSI for a given color
        bool done_with_this_colour = false;

        long s = 0;
        while (s < num_vertices && !done_with_this_colour)
        {
            long t = s+1;
            while (t < num_vertices && !done_with_this_colour)
            {
                // wanted: a (s_2, t_1) separating cut
                long s_in_D = (x_val[s][colour] > 0 && x_val[s][colour] < 1) ?
                    D_idx_of_vertex[s]+1 : D_idx_of_vertex[s];

                long t_in_D = D_idx_of_vertex[t];

                if ( instance->graph->index_matrix[s][t] < 0 &&  // non-adjacent
                     x_val[s][colour] + x_val[t][colour] > 1 &&  // might cut x*
                     s_in_D != t_in_D )                        // not contracted
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

                    if (mincut < x_val[s][colour] + x_val[t][colour] - 1 
                                 - MSI_VIOLATION_TOL)
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
                                    bool u_at_zero = x_val[u][colour] == 0;

                                    bool u_frac = x_val[u][colour] > 0 &&
                                                  x_val[u][colour] < 1 ;

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

                        violated_constr += x_vars[s][colour];
                        violated_constr += x_vars[t][colour];

                        for (vector<long>::iterator it = S.begin();
                                                    it != S.end(); ++it)
                        {
                            violated_constr += ( (-1) * x_vars[*it][colour] );
                        }

                        cuts_lhs.push_back(violated_constr);
                        cuts_rhs.push_back(1);

                        #ifdef DEBUG_MSI
                            double violating_lhs = 0;

                            cout << "### ADDED MSI: ";

                            cout << "x_" << s << " + x_" << t;
                            violating_lhs += x_val[s][colour];
                            violating_lhs += x_val[t][colour];

                            for (vector<long>::iterator it = S.begin();
                                                        it != S.end(); ++it)
                            {
                                cout << " - x_" << *it << "";
                                violating_lhs -= x_val[*it][colour];
                            }

                            cout << " <= 1 " << endl;
                            cout << right;
                            cout << setw(80) << "(lhs at current point "
                                 << violating_lhs << ")" << endl;
                            cout << left;
                        #endif

                        done_with_this_colour = true;

                        if (!SEARCH_ALL_COLOURS_FOR_MSI && !at_root_relaxation)
                            done = true;
                    }
                }

                ++t;
            }

            ++s;
        }

        this->msi_current_colour = (this->msi_current_colour+1) % num_subgraphs;
        ++colours_tried;
    }

    return (cuts_lhs.size() > 0);
}

void inline CKSCutGenerator::lift_to_minimal_separator(vector<long> &S,
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

void inline CKSCutGenerator::dfs_avoiding_set(vector<long> &S,
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

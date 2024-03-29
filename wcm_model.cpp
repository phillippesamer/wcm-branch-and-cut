#include "wcm_model.h"

/// algorithm setup switches

bool GRB_CUTS = true;             // gurobi (automatic) cuts on/off
bool GRB_HEURISTICS = true;       // gurobi (automatic) heuristics on/off
bool GRB_PREPROCESSING = true;    // gurobi (automatic) preprocessing on/off

const bool GRB_HEURISTICS_FOCUS = true; // extra focus on gurobi heuristics

const double EPSILON_TOL = 1e-5;

WCMModel::WCMModel(IO *instance)
{
    this->instance = instance;

    this->solution_weight = numeric_limits<double>::max();
    this->solution_dualbound = numeric_limits<double>::max();
    this->solution_vector_x = vector<bool>(instance->graph->num_edges, false);
    this->solution_vector_y = vector<bool>(instance->graph->num_vertices, false);
    this->solution_status = STATUS_UNKNOWN;
    this->solution_runtime = -1;

    this->lp_bound = this->lp_runtime = this->lp_passes = -1;

    try
    {
        this->env = new GRBEnv();
        this->model = new GRBModel(*env);

        create_variables();
        create_constraints();
        create_objective();

        this->cutgen = new WCMCutGenerator(model, x, y, instance);

        //model->write("wcm.lp");
    }
    catch(GRBException e)
    {
        cout << "Model construction error, code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

WCMModel::~WCMModel()
{
    delete cutgen;

    delete[] x;
    delete[] y;
    delete model;
    delete env;
}

void WCMModel::create_variables()
{
    char buffer[50];

    // binary vars x[e] = 1 iff edge e is in the matching
    x = new GRBVar[instance->graph->num_edges];
    for (long e = 0; e < instance->graph->num_edges; ++e)
    {
        sprintf(buffer, "x_%ld", e);
        x[e] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, buffer);
    }

    // binary vars y[u] = 1 iff vertex u is covered by the matching
    // NB! y is just an auxiliary alias in this formulation; no need to force integrality
    y = new GRBVar[instance->graph->num_vertices];
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        sprintf(buffer, "y_%ld", u);
        y[u] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, buffer);
    }

    model->update();
}

void WCMModel::create_constraints()
{
    ostringstream cname;

    // 1. DEGREE INEQUALITIES: x vars induce a matching
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        GRBLinExpr deg_ineq = 0;

        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];
            deg_ineq += x[e];
        }

        cname.str("");
        cname << "C1_DEGREE_" << u;
        model->addConstr(deg_ineq <= 1, cname.str());
    }

    // 2. LINKING CONSTRAINTS: y_u is the sum of x_uv for v neighbours of u
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        GRBLinExpr link_xy = 0;

        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];
            link_xy += x[e];
        }

        link_xy -= y[u];

        cname.str("");
        cname << "C2_LINKING_y" << u;
        model->addConstr(link_xy == 0, cname.str());
    }

    model->update();
}

void WCMModel::create_objective()
{
    GRBLinExpr objective_expression = 0;

    for (long e = 0; e < instance->graph->num_edges; ++e)
        objective_expression += (instance->graph->w[e]) * x[e];

    model->setObjective(objective_expression, GRB_MAXIMIZE);

    model->update();
}

int WCMModel::solve(bool logging)
{
    try
    {
        // solver features
        if (!GRB_CUTS)
            model->set(GRB_IntParam_Cuts, 0);

        if (!GRB_HEURISTICS)
            model->set(GRB_DoubleParam_Heuristics, 0);
        else
        {
            if (GRB_HEURISTICS_FOCUS)
            {
                model->set(GRB_IntParam_MIPFocus, 1);
                model->set(GRB_DoubleParam_Heuristics, 0.2);
            }
        }

        if (!GRB_PREPROCESSING)
        {
            model->set(GRB_IntParam_Presolve, 0);
            model->set(GRB_IntParam_PrePasses, 0);
            model->set(GRB_DoubleParam_PreSOS1BigM, 0);
            model->set(GRB_DoubleParam_PreSOS2BigM, 0);
            model->set(GRB_IntParam_PreSparsify, 0);
            model->set(GRB_IntParam_PreCrush, 1);
            model->set(GRB_IntParam_DualReductions, 0);
            model->set(GRB_IntParam_Aggregate, 0);
        }

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        // should disable presolve reductions that affect user cuts
        model->set(GRB_IntParam_PreCrush, 1);

        // must set parameter indicating presence of lazy constraints
        model->set(GRB_IntParam_LazyConstraints, 1);

        // trigger b&c separating blossom, minimal separator and indegree inequalities 
        model->setCallback(this->cutgen);
        model->optimize();

        if (logging)
        {
            cout << "Blossom inequalities added: "
                 <<  cutgen->blossom_counter << endl;

            cout << "Minimal separator inequalities added: "
                 <<  cutgen->minimal_separators_counter << endl;

            cout << "Indegree inequalities added: "
                 << cutgen->indegree_counter << endl;
        }

        return this->save_optimization_status();
    }
    catch(GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return 0;
    }
    catch(...)
    {
        cout << "Unexpected error during optimization inside WCMModel::solve()"
             << endl;
        return 0;
    }
}

///////////////////////////////////////////////////////////////////////////////

bool WCMModel::check_solution()
{
    /***
     * Depth-first search checking that vertices covered by the matching induce
     * a connected subgraph. NB! Assuming solution_vector_y and 
     * solution_vector_x are set.
     */

    // I. X AND Y VARS CORRECTLY LINKED

    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        if (solution_vector_y.at(u))
        {
            bool chosen_vertex_indeed_covered = false;

            for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
                 it != instance->graph->adj_list.at(u).end(); ++it)
            {
                long v = *it;
                long edge_idx = instance->graph->index_matrix[u][v];
                if (solution_vector_x.at(edge_idx))
                    chosen_vertex_indeed_covered = true;
            }

            if (!chosen_vertex_indeed_covered)
            {
                cout << endl << "ERROR (X~Y VAR LINK): chosen vertex not covered by a matching edge" << endl;
                return false;
            }
        }
        else
        {
            bool missing_vertex_indeed_not_covered = true;

            for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
                 it != instance->graph->adj_list.at(u).end(); ++it)
            {
                long v = *it;
                long edge_idx = instance->graph->index_matrix[u][v];
                if (solution_vector_x.at(edge_idx))
                {
                    missing_vertex_indeed_not_covered = false;

                    cout << "y[u=" << u << "] = " << setw(20) << fixed << setprecision(16) << y[u].get(GRB_DoubleAttr_X)
                         << "    (and solution_vector_y.at(u)=" << solution_vector_y.at(u) << ")" << endl;
                    cout << "y[v=" << v << "] = " << setw(20) << fixed << setprecision(16) << y[v].get(GRB_DoubleAttr_X)
                         << "    (and solution_vector_y.at(v)=" << solution_vector_y.at(v) << ")" << endl;
                    cout << "uv is edge #" << edge_idx << " in the graph" << endl;
                    cout << "x[" << edge_idx << "] = " << setw(20) << fixed << setprecision(16) << x[edge_idx].get(GRB_DoubleAttr_X)
                         << "    (and solution_vector_x.at(edge_idx)=" << solution_vector_x.at(edge_idx) << ")" << endl;
                }
            }

            if (!missing_vertex_indeed_not_covered)
            {
                cout << endl << "ERROR (X~Y VAR LINK): unchosen vertex is actually covered by a matching edge" << endl;
                return false;
            }
        }
    }

    // II. X VARS INDUCE A MATCHING

    vector<long> degree_in_solution = vector<long>(instance->graph->num_vertices, 0);
    bool still_a_matching = true;

    for (long e = 0; e < instance->graph->num_edges; ++e)
    {
        if (solution_vector_x.at(e))
        {
            long v1 = instance->graph->s.at(e);
            long v2 = instance->graph->t.at(e);
            degree_in_solution.at(v1)++;
            degree_in_solution.at(v2)++;

            if (degree_in_solution.at(v1) > 1 || degree_in_solution.at(v2) > 1)
                still_a_matching = false;
        }
    }

    if (!still_a_matching)
    {
        cout << endl << "ERROR (MATCHING): some vertex has degree above 1" << endl;
        return false;
    }


    // III. Y VARS INDUCE A CONNECTED SUBGRAPH

    long num_components = 0;
    vector<bool> seen = vector<bool>(instance->graph->num_vertices, false);

    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        // should enter only once, from the first covered vertex 
        if (solution_vector_y.at(u) && !seen.at(u))
        {
            dfs_to_tag_component(u, seen);
            ++num_components;
        }
    }

    if (num_components > 1)
    {
        cout << endl
             << "ERROR (CONNECTIVITY): subgraph induced by y is not connected ("
             << num_components << " components)" << endl;
        return false;
    }

    return true;
}

void WCMModel::dfs_to_tag_component(long u,
                                    vector<bool> &seen)
{
    // auxiliary dfs to check the connected component from vertex u

    seen.at(u) = true;

    for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
         it != instance->graph->adj_list.at(u).end(); ++it)
    {
        long v = *it;

        if (solution_vector_y.at(v) && !seen.at(v))
            dfs_to_tag_component(v, seen);
    }
}

///////////////////////////////////////////////////////////////////////////////

int WCMModel::save_optimization_status()
{
    /// Set class fields accordingly after call to optimize()

    this->solution_runtime = model->get(GRB_DoubleAttr_Runtime);

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        this->solution_status = AT_OPTIMUM;

        this->solution_weight = this->solution_dualbound 
                              = model->get(GRB_DoubleAttr_ObjVal);

        this->fill_solution_vectors();

        if (this->check_solution())
            cout << "PASSED SOLUTION CHECK" << endl << endl;
        else
            cout << endl
                 << "######################" << endl
                 << "### WRONG SOLUTION ###" << endl
                 << "######################" << endl << endl;

        // returning number of feasible solutions found (including sub-optimal)
        return model->get(GRB_IntAttr_SolCount);
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        this->solution_status = STATUS_UNKNOWN;

        this->solution_weight = numeric_limits<double>::max();
        this->solution_dualbound = numeric_limits<double>::max();
        this->solution_vector_x = vector<bool>(instance->graph->num_edges, false);
        this->solution_vector_y = vector<bool>(instance->graph->num_vertices, false);

        cout << "UNEXPECTED ERROR: model infeasible! (runtime "
             << solution_runtime << ")" << endl;

        return 0;
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
        this->solution_status = STATUS_UNKNOWN;

        this->solution_dualbound = model->get(GRB_DoubleAttr_ObjBound);

        if (model->get(GRB_IntAttr_SolCount) > 0)
        {
            this->solution_weight = model->get(GRB_DoubleAttr_ObjVal);

            this->fill_solution_vectors();

            if (this->check_solution())
                cout << "PASSED SOLUTION CHECK" << endl << endl;
            else
                cout << endl
                     << "######################" << endl
                     << "### WRONG SOLUTION ###" << endl
                     << "######################" << endl << endl;
        }
        else
            this->solution_weight = numeric_limits<double>::max();

        cout << "Time limit exceeded (" << solution_runtime << ")" << endl;
        cout << "Primal bound " << this->solution_weight 
             << ", dual bound " << this->solution_dualbound 
             << " (MIP gap " << 100*model->get(GRB_DoubleAttr_MIPGap) << "%)" 
             << endl;

        return 0;
    }
    else
    {
        this->solution_status = STATUS_UNKNOWN;
        this->solution_weight = numeric_limits<double>::max();
        this->solution_dualbound = numeric_limits<double>::max();
        this->solution_vector_x = vector<bool>(instance->graph->num_edges, false);
        this->solution_vector_y = vector<bool>(instance->graph->num_vertices, false);

        cout << "UNEXPECTED ERROR: unknown status after solve()" << endl;

        return 0;
    }
}

void WCMModel::fill_solution_vectors()
{
    // NB: gurobi vars are floating point, allowing +0 and -0
    this->solution_vector_x = vector<bool>(instance->graph->num_edges, false);
    this->solution_vector_y = vector<bool>(instance->graph->num_vertices, false);

    ostringstream solution_output;
    solution_output.str("");

    solution_output << "### Solution matching:" << endl;
    for (long e = 0; e < instance->graph->num_edges; ++e)
    {
        if (this->x[e].get(GRB_DoubleAttr_X) >= 0.5)
        {
            this->solution_vector_x.at(e) = true;

            solution_output << "  x_" << e << " = {";
            solution_output << instance->graph->s[e] << ", ";
            solution_output << instance->graph->t[e] << "}" << endl;
        }
    }

    solution_output << endl << "### Covered vertices: " ;

    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        if (this->y[u].get(GRB_DoubleAttr_X) >= 0.5)
        {
            this->solution_vector_y.at(u) = true;

            solution_output << u << " ";
        }
    }

    #ifdef DEBUG
        cout << endl << solution_output.str() << endl << endl;
    #endif
}

bool WCMModel::solve_lp_relax(bool logging, double time_limit, bool grb_cuts_off)
{
    /***
     * Determine the LP relaxation bound of the full IP formulation, including
     * blossom, indegree and minimal separator inequalities. Returns true iff 
     * the bound was computed successfully.
     */

    try
    {
        if (grb_cuts_off)
            model->set(GRB_IntParam_Cuts, 0);

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        // make vars continuous
        for (long e = 0; e < instance->graph->num_edges; ++e)
            x[e].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        // TO DO: no need to tell the solver that y needs to be integral, right?
        // for (long u = 0; u < instance->graph->num_vertices; ++u)
        //     y[u].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        // calculating wall clock time to solve LP relaxation
        struct timeval *clock_start
            = (struct timeval *) malloc(sizeof(struct timeval));
        struct timeval *clock_stop
            = (struct timeval *) malloc(sizeof(struct timeval));

        gettimeofday(clock_start, 0);

        // solve LP to optimality; then iterate separation and reoptimization
        model->optimize();
        bool updated = true;
        this->lp_passes = 1;
        this->lp_runtime = model->get(GRB_DoubleAttr_Runtime);

        while (updated && model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            cout << "LP relaxation pass #" << lp_passes << " (bound = "
                 << model->get(GRB_DoubleAttr_ObjVal) << ", runtime: "
                 << this->lp_runtime << ")" << endl;

            // cut generator object used only to find violated inequalities;
            // the callback in gurobi is not run in this context
            updated = cutgen->separate_lpr();

            // reoptimize
            if (updated)
            {
                model->optimize();  // includes processing pending modifications
                this->lp_passes++;

                // time limit bookkeeping
                gettimeofday(clock_stop, 0);
                unsigned long clock_time
                    = 1.e6 * (clock_stop->tv_sec - clock_start->tv_sec) +
                             (clock_stop->tv_usec - clock_start->tv_usec);

                this->lp_runtime = ((double)clock_time / (double)1.e6);

                if (this->lp_runtime > time_limit)
                {
                    cout << endl << "[LPR] Time limit exceeded" << endl;
                    updated = false;
                }
            }
        }

        free(clock_start);
        free(clock_stop);

        // loop might have broken because no violated inequality exists
        // or because the problem is actually infeasible
        if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            this->lp_bound = model->get(GRB_DoubleAttr_ObjVal);

            #ifdef DEBUG_LPR
                cout << "LPR solution: " << endl;

                ostringstream solution_output;
                solution_output.str("");
                solution_output << "(fractional) matching edges: " << endl;

                for (long e = 0; e < instance->graph->num_edges; ++e)
                {
                    if (this->x[e].get(GRB_DoubleAttr_X) > EPSILON_TOL)
                    {
                        solution_output << "    x[" << e << "] = "
                                        << x[e].get(GRB_DoubleAttr_X) << endl;
                    }
                }

                solution_output << "(fractionally) covered vertices: " << endl;

                for (long u = 0; u < instance->graph->num_vertices; ++u)
                {
                    if (this->y[u].get(GRB_DoubleAttr_X) > EPSILON_TOL)
                    {
                        solution_output << "    y[" << u << "] = "
                                        << y[u].get(GRB_DoubleAttr_X) << endl;
                    }
                }

                cout << solution_output.str() << endl;
            #endif

            cout << "[LPR] Blossom inequalities added: "
                 <<  cutgen->blossom_counter << endl;

            cout << "[LPR] Minimal separator inequalities added: "
                 <<  cutgen->minimal_separators_counter << endl;

            cout << "[LPR] Indegree inequalities added: "
                 << cutgen->indegree_counter << endl << endl;

            long x_frac = 0;
            for (long e = 0; e < instance->graph->num_edges; ++e)
            {
                double x_e = x[e].get(GRB_DoubleAttr_X);
                if (x_e > EPSILON_TOL && x_e < 1-EPSILON_TOL)
                    ++x_frac;
            }

            long y_frac = 0;
            for (long u = 0; u < instance->graph->num_vertices; ++u)
            {
                double y_u = y[u].get(GRB_DoubleAttr_X);
                if (y_u > EPSILON_TOL && y_u < 1-EPSILON_TOL)
                    ++y_frac;
            }

            if (x_frac > 0 || y_frac > 0)
            {
                cout << "[LPR] " << x_frac << " fractional x variables" << endl;
                cout << "[LPR] " << y_frac << " fractional y variables" << endl;
            }
            else
                cout << "[LPR] integer feasible solution" << endl;

            // restore IP model
            model->set(GRB_IntParam_Cuts, -1);
            for (long e = 0; e < instance->graph->num_edges; ++e)
                x[e].set(GRB_CharAttr_VType, GRB_BINARY);

            // TO DO: no need to tell the solver that y needs to be integral, right?
            // for (long u = 0; u < instance->graph->num_vertices; ++u)
            //     y[u].set(GRB_CharAttr_VType, GRB_BINARY);

            model->update();

            return true;
        }
        else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        {
            cout << "Unexpected error: LP relaxation infeasible!" << endl;
            cout << "Model runtime: " << lp_runtime << endl;
            return false;
        }
        else
        {
            cout << "Unexpected error: solve_lp_relax() got neither optimal "
                 << "nor infeasible model" << endl;
            cout << "Model runtime: " << lp_runtime << endl;
            return false;
        }
    }
    catch(GRBException e)
    {
        cout << "Error " << e.getErrorCode()
             << " during WCMModel::solve_lp_relax(): ";
        cout << e.getMessage() << endl;
        return false;
    }
    catch (...)
    {
        cout << "Unexpected error during WCMModel::solve_lp_relax()" << endl;
        return false;
    }
}

void WCMModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}

double WCMModel::get_mip_runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

double WCMModel::get_mip_gap()
{
    return model->get(GRB_DoubleAttr_MIPGap);
}

long WCMModel::get_mip_num_nodes()
{
    return model->get(GRB_DoubleAttr_NodeCount);
}

long WCMModel::get_mip_blossom_counter()
{
    return cutgen->blossom_counter;
}

long WCMModel::get_mip_indegree_counter()
{
    return cutgen->indegree_counter;
}

long WCMModel::get_mip_msi_counter()
{
    return cutgen->minimal_separators_counter;
}

#include "cks_model.h"

/// algorithm setup switches

bool CHECK_SOLUTION = true;            // check if the k subgraphs are connected
bool ORDER_COLOURS_CONSTRAINTS = true; // reduces solution symmetry

CKSModel::CKSModel(IO *instance)
{
    this->instance = instance;

    this->solution_weight = numeric_limits<double>::max();
    this->solution_dualbound = numeric_limits<double>::max();
    this->solution_vector = vector<long>(instance->graph->num_vertices, -1);
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

        this->cutgen = new CKSCutGenerator(model, x, instance);

        //model->write("cks.lp");
    }
    catch(GRBException e)
    {
        cout << "Model construction error, code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

CKSModel::~CKSModel()
{
    delete cutgen;

    for (long u=0; u < instance->graph->num_vertices; u++)
        delete[] x[u];
    delete[] x;

    delete model;
    delete env;
}

void CKSModel::create_variables()
{
    char buffer[50];

    // binary vars x[u][c] = 1 iff vertex u is assigned to subgraph c
    x = new GRBVar*[instance->graph->num_vertices];
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        x[u] = new GRBVar[instance->num_subgraphs];
        for (long c = 0; c < instance->num_subgraphs; ++c)
        {
            sprintf(buffer, "x_%ld_%ld", u, c);
            x[u][c] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, buffer);
        }
    }

    model->update();
}

void CKSModel::create_constraints()
{
    ostringstream cname;

    // 1. GUB constraint: each vertex assigned to at most one subgraph
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        GRBLinExpr gub_ineq = 0;
        for (long c = 0; c < instance->num_subgraphs; ++c)
            gub_ineq += x[u][c];

        cname.str("");
        cname << "C1_GUB_" << u;
        model->addConstr(gub_ineq <= 1, cname.str());
    }

    if (ORDER_COLOURS_CONSTRAINTS)
    {
        /***
         * Reduces symmetry by imposing that the component induced by colour c
         * is larger than that induced by colour c+1. Also avoids solutions with
         * "empty components" between non-empty ones.
         */
        for (long c = 0; c < instance->num_subgraphs - 1; ++c)
        {
            GRBLinExpr order_ineq = 0;

            // vertices in c - vertices in c+1
            for (long u = 0; u < instance->graph->num_vertices; ++u)
            {
                order_ineq += x[u][c];
                order_ineq += (-1)*x[u][c+1];
            }

            cname.str("");
            cname << "C2_ORDER_" << c;
            model->addConstr(order_ineq >= 0, cname.str());
        }
    }

    model->update();
}

void CKSModel::create_objective()
{
    GRBLinExpr objective_expression = 0;

    if (instance->recoloring_instance == true)
    {
        for (long u = 0; u < instance->graph->num_vertices; ++u)
        {
            long c = instance->original_colouring.at(u);
            double w = instance->graph->w.at(u);
            objective_expression += (w * x[u][c]);
        }

        model->setObjective(objective_expression, GRB_MAXIMIZE);
    }
    else
    {
        for (long u = 0; u < instance->graph->num_vertices; ++u)
            for (long c = 0; c < instance->num_subgraphs; ++c)
                objective_expression += (instance->graph->w[u]) * x[u][c];

        model->setObjective(objective_expression, GRB_MAXIMIZE);
    }

    model->update();
}

int CKSModel::solve(bool logging)
{
    try
    {
        // turn off all built-in cut generators
        // model->set(GRB_IntParam_Cuts, 0);

        /*
        // turn off all preprocessing and heuristics
        model->set(GRB_IntParam_Presolve, 0);
        model->set(GRB_IntParam_PrePasses, 0);
        model->set(GRB_DoubleParam_PreSOS1BigM, 0);
        model->set(GRB_DoubleParam_PreSOS2BigM, 0);
        model->set(GRB_IntParam_PreSparsify, 0);
        model->set(GRB_IntParam_PreCrush, 1);
        model->set(GRB_IntParam_DualReductions, 0);
        model->set(GRB_IntParam_Aggregate, 0);

        model->set(GRB_DoubleParam_Heuristics, 0);
        */

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        // should disable presolve reductions that affect user cuts
        model->set(GRB_IntParam_PreCrush, 1);

        // must set parameter indicating presence of lazy constraints
        model->set(GRB_IntParam_LazyConstraints, 1);

        // trigger b&c separating minimal separator and indegree inequalities 
        model->setCallback(this->cutgen);
        model->optimize();

        if (logging)
        {
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
        cout << "Unexpected error during optimization inside CKSModel::solve()"
             << endl;
        return 0;
    }
}

///////////////////////////////////////////////////////////////////////////////

bool CKSModel::check_solution()
{
    /***
     * Depth-first search checking each colour induces a connected subgraph.
     * NB! Assumes an integer feasible solution is available in vars x
     */

    vector<bool> colour_done = vector<bool>(instance->num_subgraphs, false);
    vector<bool> vertex_done
        = vector<bool>(instance->graph->num_vertices, false);

    // retrieve assignment of each vertex in the solution; -1 means uncoloured
    vector<long> colour_map = vector<long>(instance->graph->num_vertices, -1);
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        long c = 0;
        while (c < instance->num_subgraphs && colour_map.at(u) < 0)
        {
            if (this->x[u][c].get(GRB_DoubleAttr_X) > 1-EPSILON_TOL)
                colour_map.at(u) = c;
            else
                ++c;
        }
    }

    // dfs from each colour
    for (long u = 0; u < instance->graph->num_vertices; ++u)
    {
        long u_colour = colour_map.at(u);   // -1 if uncoloured

        if (vertex_done.at(u) == false && u_colour > 0)
        {
            if (colour_done.at(u_colour))
            {
                // had already processed that colour without finding u
                return false;
            }
            else
            {
                dfs_to_tag_component(u, u_colour, colour_map, vertex_done);
                colour_done.at(u_colour) = true;
            }
        }
    }

    return true;
}

void CKSModel::dfs_to_tag_component(long u,
                                    long u_colour,
                                    vector<long> &colour_map,
                                    vector<bool> &vertex_done)
{
    // auxiliary dfs to check the connected component from vertex u

    vertex_done.at(u) = true;

    for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
         it != instance->graph->adj_list.at(u).end(); ++it)
    {
        long v = *it;
        long v_colour = colour_map.at(v);

        if (v_colour == u_colour && vertex_done.at(v) == false)
            dfs_to_tag_component(v, u_colour, colour_map, vertex_done);
    }
}

///////////////////////////////////////////////////////////////////////////////

int CKSModel::save_optimization_status()
{
    /// Set class fields accordingly after call to optmize()

    //this->model->write("cks_full_model.lp");

    this->solution_runtime = model->get(GRB_DoubleAttr_Runtime);

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        this->solution_status = AT_OPTIMUM;

        this->solution_weight = this->solution_dualbound 
                              = model->get(GRB_DoubleAttr_ObjVal);

        this->solution_vector = vector<long>(instance->graph->num_vertices, -1);

        #ifdef DEBUG
            cout << "Optimal solution: " << endl;
        #endif

        ostringstream *solution_output 
            = new ostringstream[instance->num_subgraphs];

        for (long c = 0; c < instance->num_subgraphs; ++c)
        {
            solution_output[c].str("");
            solution_output[c] << "#" << c << ": ";
        }

        for (long u = 0; u < instance->graph->num_vertices; ++u)
            for (long c = 0; c < instance->num_subgraphs; ++c)
                if (this->x[u][c].get(GRB_DoubleAttr_X) > EPSILON_TOL)
                {
                    // NB: gurobi vars are floating point, allowing +0 and -0
                    this->solution_vector.at(u) = c;
                    solution_output[c] << u << " ";
                }

        #ifdef DEBUG
            for (long c = 0; c < instance->num_subgraphs; ++c)
                cout << solution_output[c].str() << endl;
        #endif

        delete[] solution_output;

       if (CHECK_SOLUTION)
        {
            if (this->check_solution())
                cout << "Passed solution check" << endl;
            else
                cout << "WRONG SOLUTION" << endl;
        }

        // returning number of feasible solutions found (including sub-optimal)
        return model->get(GRB_IntAttr_SolCount);
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        this->solution_status = STATUS_UNKNOWN;

        this->solution_weight = numeric_limits<double>::max();
        this->solution_dualbound = numeric_limits<double>::max();
        this->solution_vector = vector<long>(instance->graph->num_vertices, -1);

        cout << "UNEXPECTED ERROR: model infeasible! (runtime "
             << solution_runtime << ")" << endl;

        return 0;
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
        this->solution_status = STATUS_UNKNOWN;

        this->solution_weight = (model->get(GRB_IntAttr_SolCount) > 0) ?
                                model->get(GRB_DoubleAttr_ObjVal) :
                                numeric_limits<double>::max();
        this->solution_dualbound = model->get(GRB_DoubleAttr_ObjBound);

        cout << "Time limit exceeded (" << solution_runtime << ")" << endl;
        cout << "Dual bound " << this->solution_dualbound 
             << ", primal bound " << this->solution_weight 
             << " (MIP gap " << 100*model->get(GRB_DoubleAttr_MIPGap) << "%)" 
             << endl;

        return 0;
    }
    else
    {
        this->solution_status = STATUS_UNKNOWN;
        this->solution_weight = numeric_limits<double>::max();
        this->solution_dualbound = numeric_limits<double>::max();
        this->solution_vector = vector<long>(instance->graph->num_vertices, -1);

        cout << "UNEXPECTED ERROR: unknown status after solve()" << endl;

        return 0;
    }
}

bool CKSModel::solve_lp_relax(bool logging)
{
    /***
     * Solves the LP relaxation of the full IP formulation for connected
     * k-subpartitions, including indegree and minimal separator inequalities.
     * Returns true iff the bound was computed successfully.
     */

    try
    {
        // turn off all gurobi cut generators
        model->set(GRB_IntParam_Cuts, 0);

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        // make vars continuous
        for (long u = 0; u < instance->graph->num_vertices; ++u)
            for (long c = 0; c < instance->num_subgraphs; ++c)
                x[u][c].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        // calculating wall clock time to solve LP relaxation
        struct timeval *clock_start
            = (struct timeval *) malloc(sizeof(struct timeval));
        struct timeval *clock_stop
            = (struct timeval *) malloc(sizeof(struct timeval));

        gettimeofday(clock_start, 0);

        // solve LP to optimality; then iterate separation and reoptimization
        model->optimize();
        bool model_updated = true;
        this->lp_passes = 1;

        while (model_updated && model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            #ifdef DEBUG_LPR
            cout << "LP relaxation pass #" << lp_passes << " (bound = "
                 << model->get(GRB_DoubleAttr_ObjVal) << ")" << endl;
            #endif

            // cut generator object used only to find violated inequalities;
            // the callback in gurobi is not run in this context
            model_updated = cutgen->separate_lpr();

            // reoptimize
            if (model_updated)
            {
                model->optimize();  // includes processing pending modifications
                this->lp_passes++;
            }
        }

        // LP relaxation time
        gettimeofday(clock_stop, 0);
        unsigned long clock_time
            = 1.e6 * (clock_stop->tv_sec - clock_start->tv_sec) +
                     (clock_stop->tv_usec - clock_start->tv_usec);
        this->lp_runtime = ((double)clock_time / (double)1.e6);
        free(clock_start);
        free(clock_stop);

        // loop might have broken because no violated inequality exists
        // or because the problem is actually infeasible
        if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            this->lp_bound = model->get(GRB_DoubleAttr_ObjVal);

            #ifdef DEBUG_LPR
                cout << "LPR solution: " << endl;

                ostringstream *solution_output 
                    = new ostringstream[instance->num_subgraphs];

                for (long c = 0; c < instance->num_subgraphs; ++c)
                {
                    solution_output[c].str("");
                    solution_output[c] << "#" << c << ": " << endl;
                }

                for (long u = 0; u < instance->graph->num_vertices; ++u)
                    for (long c = 0; c < instance->num_subgraphs; ++c)
                        if (this->x[u][c].get(GRB_DoubleAttr_X) > EPSILON_TOL)
                        {
                            solution_output[c] << "    x[" << u
                                               << "][" << c << "] = "
                                               << x[u][c].get(GRB_DoubleAttr_X)
                                               << endl;
                        }

                for (long c = 0; c < instance->num_subgraphs; ++c)
                    cout << solution_output[c].str() << endl;

                delete[] solution_output;
            #endif

            if (logging)
            {
                cout << "Minimal separator inequalities added: "
                     <<  cutgen->minimal_separators_counter << endl;

                cout << "Indegree inequalities added: "
                     << cutgen->indegree_counter << endl;
            }

            // restore IP model
            model->set(GRB_IntParam_Cuts, -1);
            for (long u = 0; u < instance->graph->num_vertices; ++u)
                for (long c = 0; c < instance->num_subgraphs; ++c)
                    x[u][c].set(GRB_CharAttr_VType, GRB_BINARY);

            model->update();

            return true;
        }
        else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        {
            cout << "LP relaxation infeasible!" << endl;
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
             << " during CKSModel::solve_lp_relax(): ";
        cout << e.getMessage() << endl;
        return false;
    }
    catch (...)
    {
        cout << "Unexpected error during CKSModel::solve_lp_relax()" << endl;
        return false;
    }
}

void CKSModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}

double CKSModel::get_mip_runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

double CKSModel::get_mip_gap()
{
    return model->get(GRB_DoubleAttr_MIPGap);
}

long CKSModel::get_mip_num_nodes()
{
    return model->get(GRB_DoubleAttr_NodeCount);
}

long CKSModel::get_mip_msi_counter()
{
    return cutgen->minimal_separators_counter;
}

long CKSModel::get_mip_indegree_counter()
{
    return cutgen->indegree_counter;
}

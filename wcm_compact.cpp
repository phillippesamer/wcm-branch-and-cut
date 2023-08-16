#include "wcm_compact.h"

/// algorithm setup switches

const bool GRB_HEURISTICS_FOCUS = false; // extra focus on gurobi heuristics

const double EPSILON_TOL = 1e-5;

CompactWCMModel::CompactWCMModel(IO *instance)
{
    this->instance = instance;
    this->num_vertices = instance->graph->num_vertices;
    this->num_edges = instance->graph->num_edges;
    this->num_arcs = 2*num_edges + num_vertices;   // including arcs from artificial flow source

    this->solution_weight = numeric_limits<double>::max();
    this->solution_dualbound = numeric_limits<double>::max();
    this->solution_status = STATUS_UNKNOWN;
    this->solution_runtime = -1;

    this->solution_vector_x = vector<bool>(num_edges, false);
    this->solution_vector_y = vector<bool>(num_vertices, false);

    this->lp_bound = this->lp_runtime = -1;

    try
    {
        this->env = new GRBEnv();
        this->model = new GRBModel(*env);

        create_variables();
        create_constraints();
        create_objective();

        //model->write("wcm_compact.lp");
    }
    catch(GRBException e)
    {
        cout << "Model construction error, code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

CompactWCMModel::~CompactWCMModel()
{
    delete[] x;
    delete[] y;
    delete[] f;
    delete model;
    delete env;
}

void CompactWCMModel::create_variables()
{
    char buffer[50];

    // binary vars x[e] = 1 iff edge e is in the matching
    x = new GRBVar[num_edges];
    for (long e = 0; e < num_edges; ++e)
    {
        sprintf(buffer, "x_%ld", e);
        x[e] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, buffer);
    }

    // binary vars y[a] = 1 iff arc a is in the matching
    // NB! choosing the original edge list as the first orientation of arcs in the flow network:
    // arc idx in [0, m-1] for arcs oriented as s->t
    // arc idx in [m, 2m-1] for arcs oriented as t->s
    y = new GRBVar[num_arcs];
    for (long a = 0; a < num_arcs; ++a)
    {
        sprintf(buffer, "y_%ld", a);
        y[a] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, buffer);
    }

    // arc flow vars f[a], with values up to num_vertices (from the artificial source)
    f = new GRBVar[num_arcs];
    for (long a = 0; a < num_arcs; ++a)
    {
        sprintf(buffer, "f_%ld", a);
        f[a] = model->addVar(0.0, num_vertices, 0.0, GRB_CONTINUOUS, buffer);
    }

    model->update();
}

void CompactWCMModel::create_constraints()
{
    ostringstream cname;

    // 1. DEGREE INEQUALITIES: X VARS INDUCE A MATCHING
    for (long u = 0; u < num_vertices; ++u)
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

    // 2. X~Y LINKING CONSTRAINTS: 1 ARC INCIDENT TO U IF COVERED BY MATCHING, 0 OTHERWISE 
    for (long u = 0; u < num_vertices; ++u)
    {
        GRBLinExpr link_xy = 0;

        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];

            long v_to_u_arc_idx = instance->graph->t.at(e) == u ? e : e+num_edges;
            link_xy += y[v_to_u_arc_idx];
        }

        // artificial arc from the source
        long s_to_u_arc_idx = (2 * num_edges) + u;
        link_xy += y[s_to_u_arc_idx];

        // from the rhs: sum of (undirected) edges including u
        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];
            link_xy -= x[e];
        }

        cname.str("");
        cname << "C2_LINK_XY_VERTEX_" << u;
        model->addConstr(link_xy == 0, cname.str());
    }

    // 3. EXACTLY ONE ARC LEAVING THE ARTIFICIAL SOURCE (IF ANY)
    GRBLinExpr one_arc_from_src_cnstr = 0;
    for (long u = 0; u < num_vertices; ++u)
    {
        long s_to_u_arc_idx = (2 * num_edges) + u;
        one_arc_from_src_cnstr += y[s_to_u_arc_idx];
    }
    cname.str("");
    cname << "C3_ONE_ARC_FROM_S";
    model->addConstr(one_arc_from_src_cnstr <= 1, cname.str());

    // 4. MAY OPEN ARC LEAVING U ONLY IF THERE EXISTS AN ARC ENTERING U
    for (long u = 0; u < num_vertices; ++u)
    {
        list<long>::iterator neighbours_of_u = instance->graph->adj_list.at(u).begin();
        while (neighbours_of_u != instance->graph->adj_list.at(u).end())
        {
            // lhs: only the arc leaving u
            long neighbour = *neighbours_of_u;
            long edge_idx = instance->graph->index_matrix[u][neighbour];
            long leaving_arc_idx = instance->graph->s.at(edge_idx) == u ? edge_idx : edge_idx+num_edges;

            GRBLinExpr leave_only_if_enter_cnstr = 0;
            leave_only_if_enter_cnstr += y[leaving_arc_idx];

            // from the rhs: arcs incident to u
            for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
                 it != instance->graph->adj_list.at(u).end(); ++it)
            {
                long v = *it;
                long e = instance->graph->index_matrix[u][v];

                long v_to_u_arc_idx = instance->graph->t.at(e) == u ? e : e+num_edges;
                leave_only_if_enter_cnstr -= y[v_to_u_arc_idx];
            }

            // from the rhs: artificial arc from the source
            long s_to_u_arc_idx = (2 * num_edges) + u;
            leave_only_if_enter_cnstr -= y[s_to_u_arc_idx];

            cname.str("");
            cname << "C4_LEAVE_ONLY_IF_ENTER_VERTEX_" << u;
            model->addConstr(leave_only_if_enter_cnstr <= 0, cname.str());

            ++neighbours_of_u;
        }
    }

    // 5. POSITIVE FLOW ONLY IF ARC OPEN
    for (long a = 0; a < num_arcs; ++a)
    {
        GRBLinExpr f_only_if_y = 0;

        f_only_if_y += f[a];
        f_only_if_y -= (num_vertices * y[a]);

        cname.str("");
        cname << "C5_FLOW_ONLY_IF_OPEN_ARC_" << a;
        model->addConstr(f_only_if_y <= 0, cname.str());
    }

    // 6. TOTAL FLOW FROM ARTIFICIAL SOURCE = #VERTICES COVERED BY THE MATCHING
    GRBLinExpr flow_counting_induced_subgraph = 0;

    // lhs: flow on all arcs from s
    for (long u = 0; u < num_vertices; ++u)
    {
        long s_to_u_arc_idx = (2 * num_edges) + u;
        flow_counting_induced_subgraph += f[s_to_u_arc_idx];
    }

    // from the rhs: all x vars (each matching edge contributes 2 vertices to induced subgraph)
    for (long e = 0; e < num_edges; ++e)
        flow_counting_induced_subgraph -= (2 * x[e]);

    cname.str("");
    cname << "C6_FLOW_FROM_SOURCE";
    model->addConstr(flow_counting_induced_subgraph == 0, cname.str());

    // 7. FLOW BALANCE ON ALL VERTICES EXCEPT THE SOURCE = 1 IF COVERED BY THE MATCHING, 0 OTHERWISE

    for (long u = 0; u < num_vertices; ++u)
    {
        GRBLinExpr flow_balance = 0;

        // lhs first sum: flow entering u
        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];

            long entering_arc_idx = instance->graph->t.at(e) == u ? e : e+num_edges;
            flow_balance += f[entering_arc_idx];
        }

        // flow from the source
        long s_to_u_arc_idx = (2 * num_edges) + u;
        flow_balance += f[s_to_u_arc_idx];

        // lhs second sum: -1 times flow leaving u
        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];

            long leaving_arc_idx = instance->graph->s.at(e) == u ? e : e+num_edges;
            flow_balance -= f[leaving_arc_idx];
        }

        // from the rhs: arcs entering u (1 if covered, 0 otherwise cf. constraints 2.)
        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];

            long entering_arc_idx = instance->graph->t.at(e) == u ? e : e+num_edges;
            flow_balance -= y[entering_arc_idx];
        }

        // arc from the source
        flow_balance -= y[s_to_u_arc_idx];

        cname.str("");
        cname << "C7_FLOW_BALANCE_ON_VERTEX_" << u;
        model->addConstr(flow_balance == 0, cname.str());
    }

    model->update();
}

void CompactWCMModel::create_objective()
{
    GRBLinExpr objective_expression = 0;

    for (long e = 0; e < num_edges; ++e)
        objective_expression += (instance->graph->w[e]) * x[e];

    model->setObjective(objective_expression, GRB_MAXIMIZE);

    model->update();
}

int CompactWCMModel::solve(bool logging)
{
    try
    {
        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        if (GRB_HEURISTICS_FOCUS)
        {
            model->set(GRB_IntParam_MIPFocus, 1);
            model->set(GRB_DoubleParam_Heuristics, 0.2);
        }

        model->optimize();

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
        cout << "Unexpected error during optimization inside CompactWCMModel::solve()"
             << endl;
        return 0;
    }
}

///////////////////////////////////////////////////////////////////////////////

bool CompactWCMModel::check_solution()
{
    /***
     * Depth-first search checking that vertices covered by the matching induce
     * a connected subgraph. NB! Assuming solution_vector_y and 
     * solution_vector_x are set.
     */

    // I. X AND Y VARS CORRECTLY LINKED

    for (long u = 0; u < num_vertices; ++u)
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

    vector<long> degree_in_solution = vector<long>(num_vertices, 0);
    bool still_a_matching = true;

    for (long e = 0; e < num_edges; ++e)
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
    vector<bool> seen = vector<bool>(num_vertices, false);

    for (long u = 0; u < num_vertices; ++u)
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

void CompactWCMModel::dfs_to_tag_component(long u,
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

int CompactWCMModel::save_optimization_status()
{
    /// Set class fields accordingly after call to optimize()

    this->solution_runtime = model->get(GRB_DoubleAttr_Runtime);

    this->solution_status = STATUS_UNKNOWN;
    this->solution_weight = numeric_limits<double>::max();
    this->solution_dualbound = numeric_limits<double>::max();
    this->solution_vector_x = vector<bool>(num_edges, false);
    this->solution_vector_y = vector<bool>(num_vertices, false);

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
        cout << "UNEXPECTED ERROR: model infeasible! (runtime "
             << solution_runtime << ")" << endl;

        return 0;
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
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

        cout << "Time limit exceeded (" << solution_runtime << ")" << endl;
        cout << "Primal bound " << this->solution_weight 
             << ", dual bound " << this->solution_dualbound 
             << " (MIP gap " << 100*model->get(GRB_DoubleAttr_MIPGap) << "%)" 
             << endl;

        return 0;
    }
    else
    {
        cout << "UNEXPECTED ERROR: unknown status after solve()" << endl;
        return 0;
    }
}

void CompactWCMModel::fill_solution_vectors()
{
    // NB: gurobi vars are floating point, allowing +0 and -0
    this->solution_vector_x = vector<bool>(num_edges, false);
    this->solution_vector_y = vector<bool>(num_vertices, false);

    ostringstream solution_output;
    solution_output.str("");

    solution_output << "### Solution matching:" << endl;
    for (long e = 0; e < num_edges; ++e)
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

    for (long u = 0; u < num_vertices; ++u)
    {
        for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
             it != instance->graph->adj_list.at(u).end(); ++it)
        {
            long v = *it;
            long e = instance->graph->index_matrix[u][v];

            long v_to_u_arc_idx = instance->graph->t.at(e) == u ? e : e+num_edges;

            if (this->y[v_to_u_arc_idx].get(GRB_DoubleAttr_X) >= 0.5)
            {
                double f = this->f[v_to_u_arc_idx].get(GRB_DoubleAttr_X);

                this->solution_vector_y.at(u) = true;
                solution_output << u << "(" << f << " flow units from "
                                << v << ")" << endl;
            }
        }
    }

    // still need to find the root vertex, reached by an arc from the source
    for (long u = 0; u < num_vertices; ++u)
    {
        long s_to_u_arc_idx = (2 * num_edges) + u;

        if (this->y[s_to_u_arc_idx].get(GRB_DoubleAttr_X) >= 0.5)
        {
            double f = this->f[s_to_u_arc_idx].get(GRB_DoubleAttr_X);

            this->solution_vector_y.at(u) = true;
            solution_output << u << "(" << f << " flow units from SOURCE)" << endl;
        }
    }

    #ifdef DEBUG
        cout << endl << solution_output.str() << endl << endl;
    #endif
}

bool CompactWCMModel::solve_lp_relax(bool logging, double time_limit)
{
    /***
     * Determine the LP relaxation bound of the IP formulation.
     * Returns true iff the lp is feasible and was solved to optimality.
     */

    try
    {
        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        // make vars continuous
        for (long e = 0; e < num_edges; ++e)
            x[e].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        for (long a = 0; a < num_arcs; ++a)
            y[a].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        model->optimize();

        this->lp_bound = model->get(GRB_DoubleAttr_ObjVal);
        this->lp_runtime = model->get(GRB_DoubleAttr_Runtime);

        cout << "LP relaxation bound = " << lp_bound
             << ", runtime: " << this->lp_runtime << endl;

        if (this->lp_runtime > time_limit)
            cout << endl << "[LPR] Time limit exceeded" << endl;

        if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            #ifdef DEBUG_LPR
                cout << "LPR solution: " << endl;

                ostringstream solution_output;
                solution_output.str("");
                solution_output << "(fractional) matching edges: " << endl;

                for (long e = 0; e < num_edges; ++e)
                {
                    if (this->x[e].get(GRB_DoubleAttr_X) > EPSILON_TOL)
                    {
                        solution_output << "    x[" << e << "] = "
                                        << x[e].get(GRB_DoubleAttr_X) << endl;
                    }
                }

                solution_output << "(fractionally) covered vertices: " << endl;

                for (long u = 0; u < num_vertices; ++u)
                {
                    for (list<long>::iterator it = instance->graph->adj_list.at(u).begin();
                         it != instance->graph->adj_list.at(u).end(); ++it)
                    {
                        long v = *it;
                        long e = instance->graph->index_matrix[u][v];

                        long v_to_u_arc_idx = instance->graph->t.at(e) == u ? e : e+num_edges;

                        if (this->y[v_to_u_arc_idx].get(GRB_DoubleAttr_X) >= EPSILON_TOL)
                        {
                            double f = this->f[v_to_u_arc_idx].get(GRB_DoubleAttr_X);
                            solution_output << "    y[" << u << "] = "
                                            << y[v_to_u_arc_idx].get(GRB_DoubleAttr_X)
                                            << "  (" << f << " flow units from "
                                            << v << ")"
                                            << endl;
                        }
                    }
                }

                // still need to check arcs from the source
                for (long u = 0; u < num_vertices; ++u)
                {
                    long s_to_u_arc_idx = (2 * num_edges) + u;

                    if (this->y[s_to_u_arc_idx].get(GRB_DoubleAttr_X) >= EPSILON_TOL)
                    {
                        double f = this->f[s_to_u_arc_idx].get(GRB_DoubleAttr_X);

                        solution_output << "    y[" << u << "] = "
                                        << y[v_to_u_arc_idx].get(GRB_DoubleAttr_X)
                                        << "  (" << f << " flow units from SOURCE)"
                                        << endl;
                    }
                }

                cout << solution_output.str() << endl;
            #endif

            long x_frac = 0;
            for (long e = 0; e < num_edges; ++e)
            {
                double x_e = x[e].get(GRB_DoubleAttr_X);
                if (x_e > EPSILON_TOL && x_e < 1-EPSILON_TOL)
                    ++x_frac;
            }

            long y_frac = 0;
            for (long a = 0; a < num_arcs; ++a)
            {
                double y_a = y[a].get(GRB_DoubleAttr_X);
                if (y_a > EPSILON_TOL && y_a < 1-EPSILON_TOL)
                    ++y_frac;
            }

            if (x_frac > 0 || y_frac > 0)
            {
                cout << "[LPR] " << x_frac << " fractional x variables" << endl;
                cout << "[LPR] " << y_frac << " fractional y variables" << endl;
            }
            else
            {
                cout << "[LPR] integer feasible solution" << endl;
                this->save_optimization_status();
            }

            // restore IP model
            for (long e = 0; e < num_edges; ++e)
                x[e].set(GRB_CharAttr_VType, GRB_BINARY);

            for (long a = 0; a < num_arcs; ++a)
                y[a].set(GRB_CharAttr_VType, GRB_BINARY);

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
             << " during CompactWCMModel::solve_lp_relax(): ";
        cout << e.getMessage() << endl;
        return false;
    }
    catch (...)
    {
        cout << "Unexpected error during CompactWCMModel::solve_lp_relax()" << endl;
        return false;
    }
}

void CompactWCMModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}

double CompactWCMModel::get_mip_runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

double CompactWCMModel::get_mip_gap()
{
    return model->get(GRB_DoubleAttr_MIPGap);
}

long CompactWCMModel::get_mip_num_nodes()
{
    return model->get(GRB_DoubleAttr_NodeCount);
}

#include "cks_model.h"

CKSModel::CKSModel(IO *instance)
{
    this->instance = instance;

    this->solution_weight = numeric_limits<double>::max();
    this->solution_dualbound = numeric_limits<double>::max();
    this->solution_vector = vector<long>(instance->graph->num_vertices, -1);
    this->solution_status = STATUS_UNKNOWN;
    this->solution_runtime = -1;

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
            x[u][c] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, buffer);
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

    model->update();
}

void CKSModel::create_objective()
{
    // TO DO: in the convex recoloring application, change to maximization and
    // color-dependent weights: w[u][c] = w(u) if c is the original colour of u;
    // w[u][c] = 0 otherwise

    GRBLinExpr objective_expression = 0;

    for (long u = 0; u < instance->graph->num_vertices; ++u)
        for (long c = 0; c < instance->num_subgraphs; ++c)
            objective_expression += (instance->graph->w[u]) * x[u][c];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();
}

int CKSModel::solve(bool logging)
{
    try
    {
        /*
        // turn off all built-in cut generators?
        model->set(GRB_IntParam_Cuts, 0);

        // turn off all preprocessing and heuristics?
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
        //model->setCallback(this->cutgen);
        model->optimize();

        // additional information on oci cuts added
        if (logging)
        {
            cout << cutgen->minimal_separators_counter
                 << " minimal separator inequalities added" << endl;

            for (map<long,long>::iterator
                 it = cutgen->minimal_separators_len.begin();
                 it != cutgen->minimal_separators_len.end(); ++it)
            {
                    cout << setw(6) << it->second << " separators of size "
                         << it->first << endl;
            }

            cout << cutgen->indegree_counter
                 << " indegree inequalities added" << endl;

            for (map<long,long>::iterator
                 it = cutgen->indegree_len.begin();
                 it != cutgen->indegree_len.end(); ++it)
            {
                    cout << setw(6) << it->second <<" indegree vectors of size "
                         << it->first << endl;
            }
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
                cout << solution_output[c] << endl;
        #endif

        delete[] solution_output;

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

double CKSModel::runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

void CKSModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}


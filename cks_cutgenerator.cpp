#include "cks_cutgenerator.h"

/// algorithm setup switches

bool SEPARATE_MSI = true;               // MSI = minimal separator inequalities
bool SEPARATE_INDEGREE = true;

bool CUTS_AT_ROOT_ONLY = false;

bool STORE_MSI_CUT_POOL = false;
bool STORE_INDEGREE_CUT_POOL = false;

// at most 14 without changing everything to long double (which gurobi ignores)
#define SEPARATION_PRECISION_IN_IP 14

///////////////////////////////////////////////////////////////////////////////

CKSCutGenerator::CKSCutGenerator(GRBModel *model, GRBVar **x_vars, IO *instance)
{
    this->model = model;
    this->x_vars = x_vars;
    this->instance = instance;

    this->num_vertices = instance->graph->num_vertices;
    this->num_subgraphs = instance->num_subgraphs;

    this->minimal_separators_counter = 0;
    this->indegree_counter = 0;
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

bool CKSCutGenerator::run_minimal_separators_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    // run separation algorithm from "Partitioning a graph into balanced
    // connected classes - Formulations, separation and experiments", 2021,
    // by [Miyazawa, Moura, Ota, Wakabayashi]
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

    // TO DO: ALL

    return false;
}

bool CKSCutGenerator::run_indegree_separation(int kind_of_cut)
{
    /// wrapper for the separation procedure to suit different execution contexts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    // run separation algorithm from "On imposing connectivity constraints in
    // integer programs", 2017, by [Wang, Buchanan, Butenko]
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
    /// Solve the separation problem for indegree inequalities
















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

#ifndef _CKS_CUT_GEN_H_
#define _CKS_CUT_GEN_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

#include "gurobi_c++.h"

// using the preflow-push algorithm in COIN-OR:LEMON (see: www.lemon.cs.elte.hu)
#include <lemon/concepts/digraph.h>
#include <lemon/smart_graph.h>
#include <lemon/preflow.h>
using namespace lemon;

#include "io.h"
#include "cks_model.h"

// kinds of cuts
#define ADD_USER_CUTS 1
#define ADD_LAZY_CNTRS 2
#define ADD_STD_CNTRS 3

/***
 * \file cks_cutgenerator.h
 * 
 * Module for the Gurobi callback class, implementing the separation procedure
 * for minimal separators and indegree inequalities
 * 
 * Extends (and implements) the abstract base class in Gurobi.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 22.12.2022
 */

class CKSCutGenerator: public GRBCallback
{
public:
    CKSCutGenerator(GRBModel *, GRBVar**, IO*);
    virtual ~CKSCutGenerator() { }

protected:
    friend class CKSModel;

    void callback();
    bool separate_lpr();

    // input and model data
    IO *instance;
    GRBModel *model;
    long num_vertices;
    long num_subgraphs;
    bool at_root_relaxation;

    GRBVar **x_vars;
    double **x_val;
    void inline clean_x_val_beyond_precision(int);

    long indegree_counter;
    bool run_indegree_separation(int);
    bool separate_indegree(vector<GRBLinExpr> &, vector<long> &);

    long minimal_separators_counter;
    bool run_minimal_separators_separation(int);
    bool separate_minimal_separators(vector<GRBLinExpr> &, vector<long> &);
    void inline lift_to_minimal_separator(vector<long> &,
                                          vector<bool> &,
                                          long,
                                          long);
    void inline dfs_avoiding_set(vector<long> &,
                                 vector<bool> &,
                                 long,
                                 vector<bool> &,
                                 long &);

    // only used with SEARCH_ALL_COLOURS_FOR_MSI = false
    long msi_current_colour;
};

#endif

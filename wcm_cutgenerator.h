#ifndef _WCM_CUT_GEN_H_
#define _WCM_CUT_GEN_H_

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
#include "wcm_model.h"

// kinds of cuts
#define ADD_USER_CUTS 1
#define ADD_LAZY_CNTRS 2
#define ADD_STD_CNTRS 3

/***
 * \file wcm_cutgenerator.h
 * 
 * Module for the Gurobi callback class, implementing the separation procedure
 * for blossom, minimal separators and indegree inequalities
 * 
 * Extends (and implements) the abstract base class in Gurobi.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 02.07.2023
 */

class WCMCutGenerator: public GRBCallback
{
public:
    WCMCutGenerator(GRBModel*, GRBVar*, IO*);
    virtual ~WCMCutGenerator() { }

protected:
    friend class WCMModel;

    void callback();
    bool separate_lpr();

    // input and model data
    IO *instance;
    GRBModel *model;
    long num_vertices;
    bool at_root_relaxation;

    GRBVar *y_vars;
    double *y_val;
    void inline clean_vars_beyond_precision(int);

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
};

#endif

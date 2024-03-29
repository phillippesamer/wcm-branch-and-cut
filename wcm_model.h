#ifndef _WCM_MODEL_H_
#define _WCM_MODEL_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <sys/time.h>
#include <utility>

#include "gurobi_c++.h"

#include "io.h"
#include "wcm_cutgenerator.h"

/***
 * \file wcm_model.h
 * 
 * Module for the integer programming model to find connected matchings of
 * maximum weight using the Gurobi solver API.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 02.07.2023
 */

class WCMCutGenerator;

class WCMModel
{
public:
    WCMModel(IO*);
    virtual ~WCMModel();

    int solve(bool);
    double solution_weight;
    double solution_dualbound;
    vector<bool> solution_vector_x;   // edge in matching or not
    vector<bool> solution_vector_y;   // vertex covered by the matching or not
    ModelStatus solution_status;
    double solution_runtime;

    bool solve_lp_relax(bool, double, bool);
    double lp_bound;
    double lp_runtime;
    long lp_passes;

    void set_time_limit(double);

    // further info methods
    double get_mip_runtime();
    double get_mip_gap();
    long get_mip_num_nodes();
    long get_mip_blossom_counter();
    long get_mip_indegree_counter();
    long get_mip_msi_counter();

protected:
    IO *instance;

    GRBEnv *env;
    GRBModel *model;
    GRBVar *x;
    GRBVar *y;

    void create_variables();
    void create_constraints();
    void create_objective();

    WCMCutGenerator *cutgen;

    int save_optimization_status();
    void fill_solution_vectors();

    bool check_solution();
    void dfs_to_tag_component(long, vector<bool>&);
};

#endif

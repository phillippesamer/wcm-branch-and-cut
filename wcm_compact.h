#ifndef _WCM_COMPACT_H_
#define _WCM_COMPACT_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <sys/time.h>
#include <utility>

#include "gurobi_c++.h"

#include "io.h"

/***
 * \file wcm_compact.h
 * 
 * Module for the compact integer programming model using an arc-flow
 * formulation to find connected matchings of maximum weight using the Gurobi
 * solver API.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 04.08.2023
 */

class CompactWCMModel
{
public:
    CompactWCMModel(IO*);
    virtual ~CompactWCMModel();

    int solve(bool);
    double solution_weight;
    double solution_dualbound;
    vector<bool> solution_vector_x;   // edge in matching or not
    vector<bool> solution_vector_y;   // vertices in induced subgraph
    ModelStatus solution_status;
    double solution_runtime;

    bool solve_lp_relax(bool, double);
    double lp_bound;
    double lp_runtime;

    void set_time_limit(double);

    // further info methods
    double get_mip_runtime();
    double get_mip_gap();
    long get_mip_num_nodes();

protected:
    IO *instance;

    GRBEnv *env;
    GRBModel *model;
    GRBVar *x;
    GRBVar *y;
    GRBVar *f;

    long num_vertices;
    long num_edges;
    long num_arcs;

    void create_variables();
    void create_constraints();
    void create_objective();

    int save_optimization_status();
    void fill_solution_vectors();

    bool check_solution();
    void dfs_to_tag_component(long, vector<bool>&);
};

#endif

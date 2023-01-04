#ifndef _CKSModel_H_
#define _CKSModel_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "gurobi_c++.h"

#include "io.h"
#include "cks_cutgenerator.h"

#define EPSILON_TOL 0.00000001

enum ModelStatus {AT_OPTIMUM, STATUS_UNKNOWN};

/***
 * \file cks_model.h
 * 
 * Module for the integer programming model to find connected k subpartitions
 * of minimum weight, using the Gurobi solver API.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 22.12.2022
 */

class CKSCutGenerator;

class CKSModel
{
public:
    CKSModel(IO*);
    virtual ~CKSModel();
    
    int solve(bool);
    double solution_weight;
    double solution_dualbound;
    vector<bool> solution_vector;
    ModelStatus solution_status;
    double solution_runtime;

    double runtime();
    void set_time_limit(double);

protected:
    IO *instance;

    GRBEnv *env;
    GRBModel *model;
    GRBVar *x;

    void create_variables();
    void create_constraints();
    void create_objective();

    CKSCutGenerator *cutgen;

    int save_optimization_status();
};

#endif

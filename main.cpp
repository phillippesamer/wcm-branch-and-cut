/***
 * \file main.cpp
 * 
 * Branch-and-cut algorithm to find a maximum weight connected k-subpartition
 * in a graph (disjoint vertex subsets inducing k connected subgraphs).
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 22.12.2022
 */

#include "io.h"
#include "cks_model.h"

#include <cstdlib>
#include <fstream>

using namespace std;

// execution switches
double RUN_CKS_WITH_TIME_LIMIT = 1800;

int main(int argc, char **argv)
{
    // 0. PARSE INPUT FILE

    IO* instance = new IO();

    if (argc == 3)
    {
        instance->num_subgraphs = atol(argv[2]);

        if (instance->parse_single_weight_input_file(string(argv[1])) == false)
        {
            cout << "unable to parse input file" << endl;
            delete instance;
            return 0;
        }
    }
    else if (argc == 2)
    {
        cout << "convex recoloring" << endl;
        if (instance->parse_CR_input_file(string(argv[1])) == false)
        {
            cout << "unable to parse CR input file" << endl;
            delete instance;
            return 0;
        }
    }
    else
    {
        cout << "usage: \t" << argv[0]
             << " [input file path] [number of subgraphs]" << endl << endl;
        cout << "assuming a CR instance if the number of subgraphs is not given"
             << endl;

        delete instance;
        return 0;
    }

    // 1. BUILD AND SOLVE CORRESPONDING IP MODEL

    CKSModel *model = new CKSModel(instance);
    model->set_time_limit(RUN_CKS_WITH_TIME_LIMIT);

    //model->solve_lp_relax(true);

    model->solve(true);

    cout << "_____________________________________________________________________________" << endl << endl;

    cout << "cks weight: " << model->solution_dualbound
         << " (runtime " << fixed << model->solution_runtime << ")";
    if (model->solution_status != AT_OPTIMUM)
        cout << " *** NOT OPTIMAL ***";
    cout << endl;

    cout << "_____________________________________________________________________________" << endl << endl;
    /*
    */
    
    delete model;
    delete instance;
    return 0;
}

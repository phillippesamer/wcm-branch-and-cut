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
bool CONVEX_RECOLORING_INSTANCE = false;   // e.g. as of ITOR'2020

double RUN_CKS_WITH_TIME_LIMIT = 1800;

bool WRITE_LATEX_TABLE_ROW = true;
string LATEX_TABLE_FILE_PATH = string("xp4table.dat");

int main(int argc, char **argv)
{
    // 0. PARSE INPUT FILE

    IO* instance = new IO();

    if (argc == 3)
    {
        if (CONVEX_RECOLORING_INSTANCE)
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
            instance->num_subgraphs = atol(argv[2]);

            if ( !instance->parse_single_weight_input_file(string(argv[1])) )
            {
                cout << "unable to parse input file" << endl;
                delete instance;
                return 0;
            }
        }
    }
    else
    {
        cout << endl << "usage: \t" << argv[0]
             << " [CR instance file] [CR solution file]" << endl << endl;
        cout << "or \t" << argv[0]
             << " [input file path] [number of subgraphs]" << endl << endl;

        delete instance;
        return 0;
    }

    if (WRITE_LATEX_TABLE_ROW)
    {
        instance->save_instance_info();
        instance->save_literature_info(string(argv[2]));
    }

    // 1. BUILD AND SOLVE THE INTEGER PROGRAM

    CKSModel *model = new CKSModel(instance);

    /*
    model->solve_lp_relax(false);
    if (WRITE_LATEX_TABLE_ROW)
        instance->save_lpr_info(model->lp_bound, model->lp_runtime);
    */
    
    model->set_time_limit(RUN_CKS_WITH_TIME_LIMIT);
    model->solve(true);

    if (WRITE_LATEX_TABLE_ROW)
    {
        instance->save_ip_info(model->solution_weight,
                               model->solution_dualbound,
                               model->get_mip_gap(),
                               model->get_mip_runtime(),
                               model->get_mip_num_nodes(),
                               model->get_mip_msi_counter(),
                               model->get_mip_indegree_counter());

        instance->write_summary_info(LATEX_TABLE_FILE_PATH);
    }

    delete model;
    delete instance;
    return 0;
}

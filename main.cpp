/***
 * \file main.cpp
 * 
 * Branch-and-cut algorithm to find a maximum weight connected matching
 * in a graph (edge subset inducing a matching, whose covered vertices induce a
 * connected graph).
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 02.07.2023
 */

#include "io.h"
#include "wcm_model.h"

#include <cstdlib>
#include <fstream>

using namespace std;

// execution switches
double RUN_CKS_WITH_TIME_LIMIT = 1800;

bool WRITE_LATEX_TABLE_ROW = true;
string LATEX_TABLE_FILE_PATH = string("xp1.dat");

int main(int argc, char **argv)
{
    // 0. PARSE INPUT FILE

    IO* instance = new IO();

    if (argc == 2)
    {
        if ( !instance->parse_input_file(string(argv[1])) )
        {
            cout << "unable to parse input file" << endl;
            delete instance;
            return 0;
        }
    }
    else
    {
        cout << endl << "usage: \t" << argv[0]
             << " [input instance path]" << endl << endl;

        delete instance;
        return 0;
    }

    if (WRITE_LATEX_TABLE_ROW)
        instance->save_instance_info();

    // 1. BUILD AND SOLVE THE INTEGER PROGRAM

    WCMModel *model = new WCMModel(instance);

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

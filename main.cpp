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
#include <string>

using namespace std;

// execution switches
double RUN_WCM_WITH_TIME_LIMIT = 1800;

bool WRITE_LATEX_TABLE_ROW = true;
string LATEX_TABLE_FILE_PATH = string("xp1.dat");

int main(int argc, char **argv)
{
    // 0. PARSE INPUT FILE

    IO* instance = new IO();

    if (argc != 2)
    {
        cout << endl << "usage: \t" << argv[0]
             << " [input instance path]" << endl << endl;

        delete instance;
        return 0;
    }
    else
    {
        string file_path = string(argv[1]);
        string file_extension = file_path.substr(file_path.find_last_of(".")+1);

        bool successful_parsing = (file_extension.compare("gcc") == 0) ?
                                  instance->parse_gcc_file(file_path) :
                                  instance->parse_stp_file(file_path);

        if (!successful_parsing)
        {
                cout << "unable to parse input file" << endl;
                delete instance;
                return 0;
        }
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
    
    model->set_time_limit(RUN_WCM_WITH_TIME_LIMIT);
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

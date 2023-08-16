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
#include "wcm_compact.h"

#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

// execution switches
bool SEPARATORS_BASED_FORMULATION = false;
double RUN_WCM_WITH_TIME_LIMIT = 3600.0;

bool WRITE_LATEX_TABLE_ROW = true;
string LATEX_TABLE_FILE_PATH = string("xp15.dat");

// switches concerning the exponential formulation only
bool DEDICATED_LP_RELAXATION = true;
double DEDICATED_LPR_TIME_LIMIT = 300;
bool DEDICATED_LPR_GRB_CUTS_OFF = false;

int main(int argc, char **argv)
{
    // 1. PARSE INPUT FILE

    IO* instance = new IO();

    if (argc < 2)
    {
        cout << endl << "usage: \t" << argv[0]
             << " input_instance_path [-e]" << endl << endl;
        cout << "[-e]: flag indicating .stp format instance WITH edge weights"
             << endl << endl;

        delete instance;
        return 0;
    }
    else
    {
        bool stp_with_edge_weights = (argc > 2);

        string file_path = string(argv[1]);
        string file_extension = file_path.substr(file_path.find_last_of(".")+1);

        bool successful_parsing = (file_extension.compare("gcc") == 0) ?
                                  instance->parse_gcc_file(file_path) :
                                  instance->parse_stp_file(file_path, stp_with_edge_weights);

        if (!successful_parsing)
        {
                cout << "unable to parse input file" << endl;
                delete instance;
                return 0;
        }
    }

    /*
    if (instance->only_nonpositive_weights)
    {
        cout << endl << "*** All edges in instance " << instance->instance_id
             << " have non-positive weight." << endl;
        cout << endl << "*** Empty solution is optimal." << endl;

        return 0;
    }
    */

    if (WRITE_LATEX_TABLE_ROW)
        instance->save_instance_info();

    if (SEPARATORS_BASED_FORMULATION)
    {
        // 2.A INTEGER PROGRAM CORRESPONDING TO THE SEPARATORS-BASED FORMULATION

        WCMModel *model = new WCMModel(instance);

        if (DEDICATED_LP_RELAXATION)
        {
            model->solve_lp_relax(false,
                                  DEDICATED_LPR_TIME_LIMIT,
                                  DEDICATED_LPR_GRB_CUTS_OFF);

            if (WRITE_LATEX_TABLE_ROW)
                instance->save_lpr_info(model->lp_bound, model->lp_runtime);
        }
        
        model->set_time_limit(RUN_WCM_WITH_TIME_LIMIT - model->lp_runtime);
        model->solve(true);

        if (WRITE_LATEX_TABLE_ROW)
        {
            instance->save_bc_info(model->solution_weight,
                                   model->solution_dualbound,
                                   model->get_mip_gap(),
                                   model->get_mip_runtime(),
                                   model->get_mip_num_nodes(),
                                   model->get_mip_blossom_counter(),
                                   model->get_mip_msi_counter(),
                                   model->get_mip_indegree_counter());
            instance->write_summary_info(LATEX_TABLE_FILE_PATH);
        }

        delete model;
    }
    else
    {
        // 2.B INTEGER PROGRAM CORRESPONDING TO THE COMPACT, ARC-FLOW FORMULATION

        CompactWCMModel *model = new CompactWCMModel(instance);

        model->solve_lp_relax(false, RUN_WCM_WITH_TIME_LIMIT);
        if (WRITE_LATEX_TABLE_ROW)
            instance->save_lpr_info(model->lp_bound, model->lp_runtime);
        
        model->set_time_limit(RUN_WCM_WITH_TIME_LIMIT - model->lp_runtime);
        model->solve(true);

        if (WRITE_LATEX_TABLE_ROW)
        {
            instance->save_compact_info(model->solution_weight,
                                        model->solution_dualbound,
                                        model->get_mip_gap(),
                                        model->get_mip_runtime(),
                                        model->get_mip_num_nodes());
            instance->write_summary_info(LATEX_TABLE_FILE_PATH);
        }

        delete model;
    }

    delete instance;
    return 0;
}

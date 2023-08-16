#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <algorithm>

#include "graph.h"

using namespace std;

enum ModelStatus {AT_OPTIMUM, STATUS_UNKNOWN};

/***
 * \file io.h
 * 
 * Module for input and output functionality, including a Graph object for
 * main data structures.
 * 
 * Some classes are declared friends to avoid cumbersome get/set calls.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 02.07.2023
 */
class IO
{
public:
    IO();
    virtual ~IO();

    bool parse_gcc_file(string);
    bool parse_stp_file(string, bool);

    void save_instance_info();
    void save_lpr_info(double, double);
    void save_bc_info(double, double, double, double, long, long, long, long);
    void save_compact_info(double, double, double, double, long);
    void write_summary_info(string);

    // instance data
    string instance_id;
    string instance_id_trimmed;
    bool only_nonpositive_weights;

private:
    friend class CompactWCMModel;
    friend class WCMModel;
    friend class WCMCutGenerator;

    stringstream summary_info;  // latex table row summary

    Graph *graph;  // different representations of the original graph
};

#endif

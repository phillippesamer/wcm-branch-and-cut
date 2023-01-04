#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>

#include "graph.h"

using namespace std;

/***
 * \file io.h
 * 
 * Module for input and output functionality, including a Graph object for
 * main data structures.
 * 
 * Some classes are declared friends to avoid cumbersome get/set calls.
 * 
 * \author Phillippe Samer <samer@uib.no>
 * \date 22.12.2022
 */
class IO
{
public:
    IO();
    virtual ~IO();
    
    bool parse_input_file(string);

    // instance data
    long k;
    string instance_id;
    string instance_id_trimmed;

private:
    friend class CKSModel;
    friend class CKSCutGenerator;

    Graph *graph;  // different representations of the original graph
};

#endif

#include "io.h"

IO::IO()
{
    // only for consistency (while an input parsing method is not called)
    this->graph = new Graph();
    this->num_subgraphs = 1;
    this->recoloring_instance = false;
    this->original_colouring = vector<long>();
}

IO::~IO()
{
    delete graph;
}

bool IO::parse_single_weight_input_file(string filename)
{
    /// Simple input with a weight for each vertex (subgraphs do not matter)

    long num_vertices, num_edges;

    ifstream input_fh(filename);
    
    if (input_fh.is_open())
    {
        string line;
        
        // skip comment lines
        do
        {
            getline(input_fh, line);
        }
        while(line.find("#") != string::npos);   // there is a trail

        instance_id.assign(line);

        // trimmed instance id: contents after last slash and before last dot
        size_t dot_pos = filename.find_last_of(".");
        size_t last_slash_pos = filename.find_last_of("/\\");
        instance_id_trimmed = filename.substr(last_slash_pos+1,
                                              dot_pos-1 - last_slash_pos);

        // 2 lines for number of vertices and edges
        input_fh >> num_vertices;
        input_fh >> num_edges;

        // initialize graph (own structures only; lemon object at the end)
        delete graph;
        this->graph = new Graph(num_vertices,num_edges);
        this->graph->init_index_matrix();

        // m lines for edges
        for (long line_idx = 0; line_idx < num_edges; ++line_idx)
        {
            long i, j;
            
            input_fh >> i;
            graph->s.push_back(i);

            input_fh >> j;
            graph->t.push_back(j);
            
            graph->adj_list[i].push_back(j);
            graph->adj_list[j].push_back(i);
            
            // should never happen
            if (graph->index_matrix[i][j] >= 0 ||
                graph->index_matrix[j][i] >= 0 )
            {
                cerr << "ERROR: repeated edge in input file line "
                     << line_idx << endl << endl;
                return false;
            }
            
            // store index of current edge
            graph->index_matrix[i][j] = line_idx;
            graph->index_matrix[j][i] = line_idx;
        }

        // n lines for vertex weights
        for (long line_idx = 0; line_idx < num_vertices; ++line_idx)
        {
            double w;
            input_fh >> w;
            graph->w.push_back(w);
        }

        input_fh.close();
    }
    else
    {
        cerr << "ERROR: could not open file (might not exist)." << endl;
        return false;
    }

    // lemon data structure initialization
    graph->init_lemon();

    return true;
}

bool IO::parse_CR_input_file(string filename)
{
    /// convex recoloring instance

    recoloring_instance = true;

    long num_vertices, num_edges;

    ifstream input_fh(filename);
    
    if (input_fh.is_open())
    {
        string line;
        
        // skip comment lines
        do
        {
            getline(input_fh, line);
        }
        while(line.find("#") != string::npos);   // there is a trail

        instance_id.assign(line);

        // trimmed instance id: contents after last slash and before last dot
        size_t dot_pos = filename.find_last_of(".");
        size_t last_slash_pos = filename.find_last_of("/\\");
        instance_id_trimmed = filename.substr(last_slash_pos+1,
                                              dot_pos-1 - last_slash_pos);

        // 3 lines for number of vertices, edges, and colours
        input_fh >> num_vertices;
        input_fh >> num_edges;
        input_fh >> num_subgraphs;

        // initialize graph (own structures only; lemon object at the end)
        delete graph;
        this->graph = new Graph(num_vertices,num_edges);
        this->graph->init_index_matrix();

        // m lines for edges
        for (long line_idx = 0; line_idx < num_edges; ++line_idx)
        {
            long i, j;
            
            input_fh >> i;
            graph->s.push_back(i);

            input_fh >> j;
            graph->t.push_back(j);
            
            graph->adj_list[i].push_back(j);
            graph->adj_list[j].push_back(i);
            
            // should never happen
            if (graph->index_matrix[i][j] >= 0 ||
                graph->index_matrix[j][i] >= 0 )
            {
                cerr << "ERROR: repeated edge in input file line "
                     << line_idx << endl << endl;
                return false;
            }
            
            // store index of current edge
            graph->index_matrix[i][j] = line_idx;
            graph->index_matrix[j][i] = line_idx;
        }

        // n lines describing the original coloring (format: vertex colour)
        original_colouring.assign(num_vertices, 0);
        for (long line_idx = 0; line_idx < num_vertices; ++line_idx)
        {
            long v, c;
            input_fh >> v;
            input_fh >> c;
            original_colouring.at(v) = (c-1);   // NB! c in {1, ..., k}!

            // NB! Setting all vertex weights to 1 in the current benchmark
            graph->w.push_back(1);
        }

        input_fh.close();
    }
    else
    {
        cerr << "ERROR: could not open file (might not exist)." << endl;
        return false;
    }

    // lemon data structure initialization
    graph->init_lemon();

    return true;
}

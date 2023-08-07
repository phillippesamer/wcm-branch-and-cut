#include "io.h"

IO::IO()
{
    summary_info = stringstream();

    // only for consistency (while an input parsing method is not called)
    this->graph = new Graph();
}

IO::~IO()
{
    delete graph;
}

bool IO::parse_gcc_file(string filename)
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
            double w;
            
            input_fh >> i;
            graph->s.push_back(i);

            input_fh >> j;
            graph->t.push_back(j);

            input_fh >> w;
            graph->w.push_back(w);
            
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

bool IO::parse_stp_file(string filename, bool edge_weights_given)
{
    /***
     * Steiner Tree Problem format used in DIMACS challenge benchmark instances
     * If edge_weights_given is set to true, edge weights are expected in each
     * edge line. Otherwise, we add the weights of the endpoint vertices of each
     * edge to determine its weight.
     */

    // instance id from file name: what's after last slash and before last dot
    size_t dot_pos = filename.find_last_of(".");
    size_t last_slash_pos = filename.find_last_of("/\\");
    instance_id = filename.substr(last_slash_pos+1, dot_pos-1 - last_slash_pos);
    instance_id_trimmed = instance_id;

    long num_vertices, num_edges;

    ifstream input_fh(filename);
    
    if (input_fh.is_open())
    {
        string line, word;
        
        // 1. SKIP IDENTIFICATION LINE AND COMMENT SECTION
        do
        {
            getline(input_fh, line);
        }
        while(line.find("SECTION Graph") == string::npos);

        // 2. TWO LINES FOR NUMBER OF VERTICES AND EDGES
        input_fh >> word;          // e.g. "Nodes 2853"
        input_fh >> num_vertices;
        input_fh >> word;          // e.g. "Edges 3335"
        input_fh >> num_edges;

        // initialize graph (own structures only; lemon object at the end)
        delete graph;
        this->graph = new Graph(num_vertices,num_edges);
        this->graph->init_index_matrix();

        // 3. m LINES FOR EDGES
        for (long line_idx = 0; line_idx < num_edges; ++line_idx)
        {
            // must skip a word (the edge line marker "E"): e.g. "E 1 2" 
            input_fh >> word;

            // NB! in this format, the vertices are labeled in [1, n]
            long i, j;
            input_fh >> i;
            input_fh >> j;
            --i;
            --j;

            if (edge_weights_given)
            {
                double weight;
                input_fh >> weight;
                graph->w.push_back(weight);
            }

            graph->s.push_back(i);
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

        if (!edge_weights_given)
        {
            /**
             * 4. ADVANCE TO TERMINALS SECTION: USING VERTEX WEIGHT TO DETERMINE 
             * EDGE WEIGHTS. EACH EDGE IS ASSIGNED A WEIGHT EQUAL TO THE SUM OF THE
             * WEIGHTS OF ITS TERMINALS.
             */
            do
            {
                getline(input_fh, line);
            }
            while(line.find("SECTION Terminals") == string::npos);

            // skip line for number of terminals e.g. "Terminals 2853"
            getline(input_fh, line);

            // n lines for vertex weights (not given in order...!)
            vector<double> tmp_vertex_weights = vector<double>(num_vertices, 0.0);
            for (long line_idx = 0; line_idx < num_vertices; ++line_idx)
            {
                // must skip a word (the terminal line marker "T"): e.g. "T 2064 -10.5885201522015"
                input_fh >> word;

                // NB! Remember: in this format, the vertices are labeled in [1, n]
                long u;
                input_fh >> u;
                --u;

                double w;
                input_fh >> w;
                tmp_vertex_weights.at(u) = w;
            }

            for (long e = 0; e < num_edges; ++e)
            {
                long v1 = graph->s.at(e);
                long v2 = graph->t.at(e);
                
                double w = tmp_vertex_weights.at(v1) + tmp_vertex_weights.at(v2);
                
                graph->w.push_back(w);
            }
        }

        #ifdef DEBUG
            cout << endl << "### STP instance parsed";
            cout << endl << "  id = " << instance_id_trimmed
                 << endl << "  n  = " << num_vertices
                 << endl << "  m  = " << num_edges
                 << endl << "  edge 0  = {" << graph->s.at(0) << ", " << graph->t.at(0) << "}"
                 << endl << "  edge " << graph->s.size() << "  = {" << graph->s.back() << ", " << graph->t.back() << "}"
                 //<< endl << "  tmp_weight of v0 = " << tmp_vertex_weights.at(0)
                 //<< endl << "  tmp_weight of v" << num_vertices-1 << " = " << tmp_vertex_weights.at(num_vertices-1)
                 << endl << "  weight of edge 0 = " << graph->w.at(0)
                 << endl << "  weight of edge " << graph->w.size() << " = " << graph->w.back()
                 << endl << endl;
        #endif

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

void IO::save_instance_info()
{
    /// save instance info: id  n  m  k
    summary_info << left;
    summary_info << setw(50) << instance_id_trimmed;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << graph->num_vertices;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << graph->num_edges;
    summary_info << setw(8) << "  &&  ";
}

void IO::save_lpr_info(double lp_bound, double lp_time)
{
    /// save lp relaxation info: bound time
    summary_info << setw(8) << fixed << setprecision(2) << lp_bound;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << lp_time;
    summary_info << setw(8) << "  &&  ";
}

void IO::save_bc_info(double lb,
                      double ub,
                      double gap,
                      double time,
                      long node_count,
                      long blossom_count,
                      long msi_count,
                      long indegree_count)
{
    /// save mip info: lb ub gap time #nodes #blossom #msi #indegree

    summary_info << setw(8) << fixed << setprecision(2) << lb;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << ub;

    summary_info << setw(8) << "  &  ";
    if (gap < 10) // < 1000%
    {
        double percentual_gap = 100 * gap;
        summary_info << setw(8) << fixed << setprecision(2) << percentual_gap;
    }
    else
        summary_info << setw(8) << " -- ";

    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << time;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << node_count;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << blossom_count;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << msi_count;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << indegree_count;
    summary_info << setw(8) << "  \\\\  ";
}

void IO::save_compact_info(double lb,
                           double ub,
                           double gap,
                           double time,
                           long node_count)
{
    /// save mip info: lb ub gap time #nodes

    summary_info << setw(8) << fixed << setprecision(2) << lb;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << ub;

    summary_info << setw(8) << "  &  ";
    if (gap < 10) // < 1000%
    {
        double percentual_gap = 100 * gap;
        summary_info << setw(8) << fixed << setprecision(2) << percentual_gap;
    }
    else
        summary_info << setw(8) << " -- ";

    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << time;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << node_count;
    summary_info << setw(8) << "  \\\\  ";
}

void IO::write_summary_info(string output_file_path)
{
    /// write all the saved info as a line in the given file

    ofstream xpfile(output_file_path.c_str(), ofstream::app);
    if (xpfile.is_open())
    {
        xpfile << summary_info.str();
        xpfile << endl;
        xpfile.close();
    }
    else
    {
        cout << "ERROR: unable to write XP file; dumping to screen:" << endl;
        cout << summary_info.str();
    }
}

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

bool IO::parse_input_file(string filename)
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

        // TO DO: edge weights instead
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

void IO::save_instance_info()
{
    /// save instance info: id  n  m  k
    summary_info << left;
    summary_info << setw(15) << instance_id_trimmed;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << graph->num_vertices;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << graph->num_edges;
    summary_info << setw(8) << "  &&  ";

    #ifdef DEBUG
        cout << "save_instance_info got: " << endl;
        cout << summary_info.str() << endl;
    #endif
}

void IO::save_lpr_info(double lp_bound, double lp_time)
{
    /// save lp relaxation info: bound time
    summary_info << setw(8) << fixed << setprecision(2) << lp_bound;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << lp_time;
    summary_info << setw(8) << "  &&  ";

    #ifdef DEBUG
        cout << "save_lpr_info got: " << endl;
        cout << summary_info.str() << endl;
    #endif
}

void IO::save_ip_info(long lb,
                      long ub,
                      double gap,
                      double time,
                      long node_count,
                      long msi_count,
                      long indegree_count)
{
    /// save mip info: lb ub gap time #nodes #msi #indegree

    summary_info << setw(8) << lb;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << ub;

    double percentual_gap = 100 * gap;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << percentual_gap;

    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << fixed << setprecision(2) << time;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << node_count;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << msi_count;
    summary_info << setw(8) << "  &  ";
    summary_info << setw(8) << indegree_count;
    summary_info << setw(8) << "  \\\\  ";

    #ifdef DEBUG
        cout << "save_ip_info got: " << endl;
        cout << summary_info.str() << endl;
    #endif
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

//     Copyright (c) 2006-2010 Hartmut Kaiser. Distributed under the Boost
//     Software License, Version 1.0. (See accompanying file
//     LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*  Class to access various ADCIRC output files "fort.63/64"
                Modified: 2 August, 2006
*/

#ifndef OUTPUT_FILE
#define OUTPUT_FILE

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#include "exception.hpp"

/****************************************
    Class to read in output files
****************************************/
namespace adcirc
{
    class output_file
    {
        // CANNOT READ AND PARSE FIRST LINE in output files .. 
        // cannot see how the format is ...
//         std::string run_description; // RUNDES = alphanumeric run description
//         std::string run_id;          // RUNID = alphanumeric run description 2
//         std::string grid_id;         // AGRID = alpha-numeric grid identification 

        unsigned int no_of_datasets;    //NDSET = the number of data sets to be written to fort.* files
        unsigned int no_of_nodes;       // NP
        unsigned int ts;                // timestep
        int record_type;                //IRTYPE = the record type 
                                        // (= 1 for elevation files, = 2 for velocity files, and = 3 for 3D velocity file)

        double model_time;             // TIME - model time (in seconds)
        unsigned int model_time_step_number;    //IT - model time step number since the beginning of the model run.

        std::vector<double> result;

        // To calculate minimum and maximum value
        double min_value, max_value;
        int min_value_node, max_value_node;

        std::string templine;

        std::string filename;
        boost::iostreams::file_source fort_in;
        boost::iostreams::filtering_istream fort_file;

        int current_step;

        std::ostream* logfile;

#if defined(USE_ADCIRC)
        // netcdf functionality
        bool is_netcdf;

        void read_netcdf_header();
        void read_netcdf(int step_min, int step_max, 
            std::vector<std::string> const& field_names, bool verbose = false);
        bool read_netcdf_next(std::vector<std::string> const& field_names, 
            std::vector<std::vector<double> >& data, 
            std::string& header, bool verbose = false);
#endif

        void skip_to(int step_min, bool verbose = false);
        void read_minmax(int step_max, int index1, int index2, 
            bool verbose = false);

    public:
        /*************************************************************
         Function to fetch the information based on various arguments
        **************************************************************/
#if defined(USE_ADCIRC)
        output_file(std::string const&, bool verbose, 
            std::vector<std::string> const& field_names, int step_min = -1, 
            int step_max = -1, std::ostream* logfile = NULL);
#endif

        output_file(std::string const&, bool verbose = false, 
            std::ostream* logfile = NULL);

        /***************************************
            To print the result values
        ****************************************/
        void print_result(void);

        /***********************************************
            To access the result (returns result)
        ***********************************************/
        std::vector<double>& get_result(void);

        /********************************************
          Access maximum and minimum value
        ********************************************/
        double get_min_value(void);
        double get_max_value(void);

        int get_min_value_node(void);
        int get_max_value_node(void);

        bool read_next(std::vector<std::string> const& field_names, 
            std::vector<std::vector<double> >& data, 
            std::string& header, bool verbose = false);

        int get_num_datasets() const { return no_of_datasets; }
        int get_timestep() const { return ts; }
    };

} // end of namespace ADCIRC

#endif

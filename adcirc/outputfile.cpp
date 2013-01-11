//     Copyright (c) 2006-2010 Hartmut Kaiser. Distributed under the Boost
//     Software License, Version 1.0. (See accompanying file
//     LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Class to access various ADCIRC output files "fort.63/64"  (function definitions)

#include <float.h>
#include <cmath>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem/convenience.hpp>

#if defined(USE_ADCIRC)
#include <adcirc/netcdf.hpp>
#endif
#include "outputfile.hpp"

inline double sqr(double x)
{
    return x * x;
}

///////////////////////////////////////////////////////////////////////////////
//            "ADCIRC_output_file" class

void adcirc::output_file::skip_to(int step_min, bool verbose)
{
    if (!fort_in.is_open())
        throw adcirc::exception("Unable to open the file: " + filename);

    for (int i = 0; i < step_min-1; ++i)
    {
        if (verbose && logfile) {
            (*logfile) << "Skipping timestep (" << i << ")";
            if (!templine.empty())
                (*logfile) << ": " << templine;
            (*logfile) << std::endl;
        }

        unsigned int no_of_timestep_nodes = no_of_nodes;

        // figure out the type of the header line
        {    
            std::stringstream strstrm(templine);
            double time_since_cold_start;
            int internal_timestep;
            double default_value = -99999.0;

            strstrm >> time_since_cold_start >> internal_timestep;
            strstrm >> no_of_timestep_nodes >> default_value;
        }

        // bypass all the values that I don't want till we reach our requested timestep
        for(unsigned int j = 0; j < no_of_timestep_nodes; ++j)
        {
            // there is one more line + no_of_nodes for every timestep
            std::string line;
            std::getline(fort_file, line);
        }

        // read first line of next block
        std::getline(fort_file, templine);
    }

    current_step = step_min;
}

#if defined(USE_ADCIRC)
void adcirc::output_file::read_netcdf_header()
{
    // handle netcdf
    try {
        // open netCDF objects
        netCDF::NcFile ncfile(filename, netCDF::NcFile::read);

        // read dimensions
        no_of_nodes = (boost::uint32_t)read_file_dimension(ncfile, "node");
        no_of_datasets = (boost::uint32_t)read_file_dimension(ncfile, "time");
        record_type = 3;
        
        is_netcdf = true;
    }
    catch (netCDF::exceptions::NcException const& e) {
        std::string msg("Caught NetCDF exception: ");
        msg += e.what();
        msg += " (file: " + filename + ")"; 
        if (logfile)
            (*logfile) << msg << std::endl;
        throw adcirc::exception(msg);
    }
}

void adcirc::output_file::read_netcdf(
    int step_min, int step_max, std::vector<std::string> const& field_names, 
    bool verbose)
{
    // handle netcdf
    try {
        // open netCDF objects
        netCDF::NcFile ncfile(filename, netCDF::NcFile::read);

        max_value = -DBL_MAX;
        min_value = DBL_MAX;
        max_value_node = -1;
        min_value_node = -1;

        result.resize(no_of_nodes);
        std::fill(result.begin(), result.end(), -99999);

        if (-1 != step_min)
            current_step = step_min-1;
        if (-1 == step_max)
            step_max = no_of_datasets;

        for (int steps = current_step; steps < step_max; ++steps/*, ++current_step*/) {
            std::vector<std::vector<double> > data;
            std::string header;
            if (!read_netcdf_next(field_names, data, header, verbose)) {
                throw adcirc::exception("Failed to read timestep (" +
                    boost::lexical_cast<std::string>(steps) + ")");
            }

            // reconcile this time step into the result
            for (std::size_t idx = 0; idx < no_of_nodes; ++idx) {
                double value = -99999;
                if (field_names.size() == 1)
                    value = data[0][idx];
                else if (data[0][idx] != -99999 && data[1][idx] != -99999)
                    value = std::sqrt(sqr(data[0][idx]) + sqr(data[1][idx]));

                if (value != -99999) {
                    if (result[idx] == -99999)
                        result[idx] = value;
                    else
                        result[idx] = (std::max)(result[idx], value);   // collect maximum

                    // for calculating minimum & maximum value
                    if (value < min_value) {
                        min_value = value;
                        min_value_node = idx;
                    }
                    if (value > max_value) {
                        max_value = value;
                        max_value_node = idx;
                    }
                }
            }
        }
    }
    catch (netCDF::exceptions::NcException const& e) {
        std::string msg("Caught NetCDF exception: ");
        msg += e.what();
        msg += " (file: " + filename + ")"; 
        if (logfile)
            (*logfile) << msg << std::endl;
        throw adcirc::exception(msg);
    }
}

/*************************************************************
 Function to fetch the information based on various arguments
**************************************************************/
// index - which column to show; 
// step_number - which timestep to choose from the whole file
adcirc::output_file::output_file(std::string const& filename, bool verbose, 
      std::vector<std::string> const& field_names, int step_min, int step_max, 
      std::ostream* logfile) 
  : min_value(DBL_MAX), max_value(-DBL_MAX), 
    min_value_node(-1), max_value_node(-1), 
    record_type(0), 
    filename(filename), fort_in(filename), current_step(0), logfile(logfile),
    is_netcdf(false)
{
    std::string ext(boost::filesystem::extension(filename));
    if (ext == ".nc" || ext == ".cdf") {
        fort_in.close();
        read_netcdf_header();
        read_netcdf(step_min, step_max, field_names, verbose);
        return;
    }

    // handle plain or compressed ascii 
    if (ext == ".gz") {
        fort_in.close();
        fort_in.open(filename, std::ios_base::in | std::ios_base::binary);
        fort_file.push(boost::iostreams::gzip_decompressor());
    }
    if(!fort_in.is_open())
        throw adcirc::exception("Unable to open the file: " + filename);

    int index1 = -1;
    int index2 = -1;
    
    if (!field_names.empty())
        index1 = boost::lexical_cast<int>(field_names[0]);
    if (field_names.size() > 1)
        index2 = boost::lexical_cast<int>(field_names[1]);

    fort_file.push(fort_in);

    if (verbose && logfile)
        (*logfile) << "Reading ADCIRC output file: " << filename << std::endl; 

    std::getline(fort_file, templine, '\n');
    std::getline(fort_file, templine, '\n');

    record_type = -1;
    if (templine.find_first_not_of(" ") != std::string::npos) {
    // this is not the maxsurge.63 file
        double temp1;

        std::stringstream strstrm(templine);
        strstrm >> no_of_datasets >> no_of_nodes >> ts >> temp1 >> record_type;
        std::getline(fort_file, templine, '\n');   // fort.63
    }
    else {
    // maxsurge.63 has one empty leading line 
        fort_file >> no_of_datasets >> no_of_nodes; 
        std::getline(fort_file, templine, '\n');     // skip that line
        record_type = 3;
    }

    if (verbose && logfile)
        (*logfile) << "Number of time steps in file: " << no_of_datasets << std::endl; 

    // we need all datasets including the last one
    if (-1 == step_max)
        step_max = no_of_datasets-1;

    if (!no_of_nodes)
        throw adcirc::exception("ZERO read no. of nodes");

    /* Go to the timestep needed and then load the values in "index" column in a std::vector */        
    // check if that timestep is valid
    if(-1 != step_min && step_min > (int)no_of_datasets) {
        throw adcirc::exception("invalid (minimum) timestamp requested: " +
            boost::lexical_cast<std::string>(step_min));
    }
    if(-1 != step_max && step_max > (int)no_of_datasets) {
        throw adcirc::exception("invalid (maximum) timestamp requested: " +
            boost::lexical_cast<std::string>(step_max));
    }

    // now check if the "index" requested is valid
    if (-1 != record_type && index1 > record_type) {
        throw adcirc::exception("invalid index1 (column data index) requested: " +
            boost::lexical_cast<std::string>(index1));
    }
    if (-1 != record_type && -1 != index2 && index2 > record_type) {
        throw adcirc::exception("invalid index2 (column data index) requested: " +
            boost::lexical_cast<std::string>(index2));
    }

    // Time step asked is valid
    skip_to(step_min, verbose);

    // now we are reading the time steps we want
    read_minmax(step_max, index1, index2, verbose);
}
#endif

///////////////////////////////////////////////////////////////////////////////
adcirc::output_file::output_file(std::string const& filename, bool verbose, 
        std::ostream* logfile) 
  : min_value(DBL_MAX), max_value(-DBL_MAX), 
    min_value_node(-1), max_value_node(-1), 
    record_type(0), 
    filename(filename), fort_in(filename), current_step(0), logfile(logfile)
#if defined(USE_ADCIRC)
  , is_netcdf(false)
#endif
{
    std::string ext(boost::filesystem::extension(filename));
#if defined(USE_ADCIRC)
    if (ext == ".nc" || ext == ".cdf") {
        fort_in.close();
        read_netcdf_header();
        return;
    }
#endif

    // handle plain or compressed ascii 
    if (ext == ".gz") {
        fort_in.close();
        fort_in.open(filename, std::ios_base::in | std::ios_base::binary);
        fort_file.push(boost::iostreams::gzip_decompressor());
    }
    if(!fort_in.is_open())
        throw adcirc::exception("Unable to open the file: " + filename);

    fort_file.push(fort_in);

    if (verbose && logfile)
        (*logfile) << "Reading ADCIRC output file: " << filename << std::endl; 

    std::getline(fort_file, templine, '\n');
    std::getline(fort_file, templine, '\n');

    record_type = -1;
    if (templine.find_first_not_of(" ") != std::string::npos) {
    // this is not the maxsurge.63 file
        double temp1;

        std::stringstream strstrm(templine);
        strstrm >> no_of_datasets >> no_of_nodes >> ts >> temp1 >> record_type;
        std::getline(fort_file, templine, '\n');   // fort.63
    }
    else {
    // maxsurge.63 has one empty leading line 
        fort_file >> no_of_datasets >> no_of_nodes; 
        std::getline(fort_file, templine, '\n');     // skip that line
        record_type = 3;
    }

    if (verbose && logfile)
        (*logfile) << "Number of time steps in file: " << no_of_datasets << std::endl; 

    if (!no_of_nodes)
        throw adcirc::exception("ZERO read no. of nodes");
}

///////////////////////////////////////////////////////////////////////////////
void adcirc::output_file::read_minmax(int step_max, int index1, int index2, 
    bool verbose)
{
    if (!fort_in.is_open())
        throw adcirc::exception("Unable to open the file: " + filename);

    max_value = -DBL_MAX;
    min_value = DBL_MAX;
    max_value_node = -1;
    min_value_node = -1;

    int max_index = index1;
    if (index2 != -1 && index2 > index1)
        max_index = index2;

    result.resize(no_of_nodes);
    std::fill(result.begin(), result.end(), -99999);

    for (int steps = current_step; steps <= step_max; ++steps, ++current_step) {
        if (verbose && logfile) {
            (*logfile) << "Reading timestep (" << steps << ")";
            if (!templine.empty())
                (*logfile) << ": " << templine;
            (*logfile) << std::endl;
        }

        unsigned int no_of_timestep_nodes = no_of_nodes;
        double default_value = -99999.0;

        // figure out the type of the header line
        {    
            std::stringstream strstrm(templine);
            double time_since_cold_start;
            int internal_timestep;

            strstrm >> time_since_cold_start >> internal_timestep;
            strstrm >> no_of_timestep_nodes >> default_value;
        }

        /* Read in the requested column name referenced by "index" */
        for(unsigned int j = 0; j < no_of_timestep_nodes; j++)
        {
            std::getline(fort_file, templine);
            if (!fort_file.good()) {
                if (logfile) 
                    (*logfile) << "reached end of file, canceling ..." << std::endl;
                return;
            }
            std::stringstream strstrm(templine);

            int idx = 0;
            strstrm >> idx; // read in the node number at each line
            
            if (idx <= 0 || idx > (int)no_of_nodes) {
                if (logfile) {
                    (*logfile) << "Invalid node number: " << idx 
                               << " (line: " << j << ")" << std::endl;
                }
                return;
            }
            --idx;

            // skip various values that we don't need in the column
            double value1 = default_value, value2 = default_value;
            for(int k = 1; k <= max_index && strstrm.good(); k++)
            {
                double value = default_value;
                strstrm >> value;

                if (strstrm.fail()) 
                    break;

                if (k == index1)
                    value1 = value;
                else if (k == index2)
                    value2 = value;
            }

            double value = 0;
            if (-1 == index2)
                value = value1;
            else if (value1 != default_value && value2 != default_value)
                value = std::sqrt(sqr(value1) + sqr(value2));
            else
                value = default_value;

            if (value != default_value) {
                if (result[idx] == default_value)
                    result[idx] = value;
                else
                    result[idx] = (std::max)(result[idx], value);   // collect maximum

                // for calculating minimum & maximum value
                if (value < min_value) {
                    min_value = value;
                    min_value_node = idx;
                }
                if (value > max_value) {
                    max_value = value;
                    max_value_node = idx;
                }
            }
        }

        // skip one more line
        std::getline(fort_file, templine);
        if (!fort_file.good() && steps != step_max) {
            if (logfile) {
                (*logfile) << "Canceling ... (reached end of file on timestep "
                           << steps << ", expected " << step_max << " steps)" 
                           << std::endl;
            }
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
#if defined(USE_ADCIRC)
namespace util
{
    inline std::size_t 
    get_dimension(netCDF::NcVar const& var, std::string const& name, std::size_t& size)
    {
        std::size_t dims = var.getDimCount();
        std::size_t dimension = 0;
        for (/**/; dimension < dims; ++dimension) {
            netCDF::NcDim dim(var.getDim(dimension));
            if (dim.getName() == name) {
                size = dim.getSize();
                break;
            }
        }

        if (dimension == dims) {
            throw adcirc::exception(
                "Couldn't find dimension '" + name + "' for requested variable '" + 
                var.getName() + '"');
        }
        return dimension;
    }
}

bool adcirc::output_file::read_netcdf_next(
    std::vector<std::string> const& field_names, 
    std::vector<std::vector<double> >& data, std::string& header, 
    bool verbose)
{
    if (!is_netcdf)
        throw adcirc::exception("Input is not a netcdf file: " + filename);

    if (verbose && logfile) {
        (*logfile) << "Reading netcdf timestep (" << current_step << ")";
        if (!templine.empty())
            (*logfile) << ": " << templine;
        (*logfile) << std::endl;
    }

    data.resize(field_names.size());
    for(unsigned int k = 0; k < field_names.size(); ++k) 
        data[k].resize(no_of_nodes);

    std::string field;
    try {
        // open netCDF objects
        netCDF::NcFile ncfile(filename, netCDF::NcFile::read);

        int i = 0;
        BOOST_FOREACH(std::string const& f, field_names)
        {
            field = f;
            netCDF::NcVar var(ncfile.getVar(field));
            std::size_t dims = var.getDimCount();
            if (dims != 2 && dims != 1) {
                throw adcirc::exception(
                    "Unexpected number of dimensions for variable '" + 
                    field + std::string("' (") + 
                    boost::lexical_cast<std::string>(dims) + std::string(")"));
            }

            std::size_t max_timesteps = no_of_datasets;
            std::size_t max_nodes = 0;
            std::size_t time_dimension = 0;
            std::size_t node_dimension = util::get_dimension(var, "node", max_nodes);

            if (dims == 2)
                time_dimension = util::get_dimension(var, "time", max_timesteps);

            if (current_step < 0 || current_step >= (int)max_timesteps) {
                throw adcirc::exception(
                    "Requested not existing timestep for variable '" + 
                    field + std::string("' (") + 
                    boost::lexical_cast<std::string>(current_step) + 
                    std::string(")"));
            }

            double default_value = -99999.0;
            double fillvalue = -99999.0;
            std::map<std::string, netCDF::NcVarAtt> attrs (var.getAtts());
            typedef std::map<std::string, netCDF::NcVarAtt>::iterator iterator;
            
            iterator it = attrs.find("missing_value");
            if (it != attrs.end()) 
                (*it).second.getValues(&default_value);

            it = attrs.find("_FillValue");
            if (it != attrs.end()) 
                (*it).second.getValues(&fillvalue);

            std::fill(data[i].begin(), data[i].end(), default_value);
            read_netcdf_var(ncfile, field, dims, time_dimension, current_step, 
                node_dimension, data[i]);

            for (int j = 0; j < data[i].size(); ++j) {
                if (data[i][j] == fillvalue)
                    data[i][j] = default_value;
            }

            ++i;
        }

        ++current_step;
    }
    catch (netCDF::exceptions::NcException const& e) {
        std::string msg("Caught NetCDF exception: ");
        msg += e.what();
        if (!field.empty())
            msg += ": " + field;
        msg += " (file: " + filename + ")"; 
        if (logfile)
            (*logfile) << msg << std::endl;
        throw adcirc::exception(msg);
    }

    return true;
}
#endif

///////////////////////////////////////////////////////////////////////////////
bool adcirc::output_file::read_next(std::vector<std::string> const& field_names, 
    std::vector<std::vector<double> >& data, std::string& header, bool verbose)
{
#if defined(USE_ADCIRC)
    if (is_netcdf) 
        return read_netcdf_next(field_names, data, header, verbose);
#endif

    if (!fort_in.is_open())
        throw adcirc::exception("Unable to open the file: " + filename);

    int index1 = -1;
    int index2 = -1;
    int index3 = -1;
    int column_count = 0;

    if (!field_names.empty()) 
        index1 = boost::lexical_cast<int>(field_names[column_count++]);
    if (field_names.size() > 1)
        index2 = boost::lexical_cast<int>(field_names[column_count++]);
    if (field_names.size() > 2)
        index3 = boost::lexical_cast<int>(field_names[column_count++]);

    int columns = (std::max)(index1, (std::max)(index2, index3));
    if (field_names.size() > 3 || columns > 3)
        throw adcirc::exception("Can't read more than 3 columns from ACIRC output file.");
    if (columns == -1)
        throw adcirc::exception("Need to specify at least one column to read from ACIRC output file.");

    if (verbose && logfile) {
        (*logfile) << "Reading timestep (" << current_step << ")";
        if (!templine.empty())
            (*logfile) << ": " << templine;
        (*logfile) << std::endl;
    }

    header = templine;
    unsigned int no_of_timestep_nodes = no_of_nodes;
    double default_value = -99999.0;

    // figure out the type of the header line
    {    
        std::stringstream strstrm(templine);
        double time_since_cold_start;
        int internal_timestep;

        strstrm >> time_since_cold_start >> internal_timestep;
        strstrm >> no_of_timestep_nodes >> default_value;
    }

    data.resize(column_count);
    for (int k = 0; k < column_count; ++k) {
        data[k].resize(no_of_nodes);
        std::fill(data[k].begin(), data[k].end(), default_value);
    }

    // Read in the requested number of columns
    for(unsigned int j = 0; j < no_of_timestep_nodes; ++j)
    {
        std::getline(fort_file, templine);
        if (!fort_file.good()) {
            if (logfile) 
                (*logfile) << "reached end of file, canceling ..." << std::endl;
            return false;
        }
        std::stringstream strstrm(templine);

        int idx = 0;
        strstrm >> idx;     // read in the node number at each line
        
        if (idx <= 0 || idx > (int)no_of_nodes) {
            if (logfile) {
                (*logfile) << "Invalid node number: " << idx 
                           << " (line: " << j << ")" << std::endl;
            }
            return false;
        }
        --idx;

        // skip various values that we don't need in the column
        double values[3] = { default_value, default_value, default_value };
        int k = 0;
        for(/**/; k < columns && strstrm.good(); ++k)
        {
            double value = default_value;
            strstrm >> value;

            if (strstrm.fail()) 
                break;

            values[k] = value;
        }

        // store data in result vector
        for (int k1 = 0, i = 0; k1 < k; ++k1) {
            if (k1 == index1 || k1 == index2 || k1 == index3)
                data[i++][idx] = values[k1];
        }            
    }

    ++current_step;

    // skip one more line
    std::getline(fort_file, templine);

    return true;
}

/***************************************
   To print the result value
****************************************/
void adcirc::output_file::print_result(void)
{
    std::cout << "\nMODEL TIME : " << model_time <<"\nTIME STEP NUMBER : " << model_time_step_number << std::endl;
    std::cout << "\nstd::vector values are : (just printing first 4 values)\n";
    int i=1;
    for(std::vector<double>::iterator iter = result.begin(); i < 5/*iter !=result.end() */; iter ++,i++)
        std::cout << i << "   :    " << (*iter)<< std::endl;
}

/***********************************************
  To access the result (returns result)
***********************************************/
std::vector<double>& adcirc::output_file::get_result() 
{
    return result; 
}

/*************************************
  Access minimum & maximum value
*************************************/
double adcirc::output_file::get_min_value()
{
    return min_value;
}

double adcirc::output_file::get_max_value()
{
    return max_value;
}

int adcirc::output_file::get_min_value_node()
{
    return min_value_node;
}

int adcirc::output_file::get_max_value_node()
{
    return max_value_node;
}

// end of class "ADCIRC_output_file" functions


//     Copyright (c) 2006 Hartmut Kaiser. Distributed under the Boost
//     Software License, Version 1.0. (See accompanying file
//     LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem/convenience.hpp>

#include "gridfile.hpp"
#include "exception.hpp"
//#include "netcdf.hpp"

///////////////////////////////////////////////////////////////////////////////
namespace adcirc
{

///////////////////////////////////////////////////////////////////////////////
// Constructor function 
grid_file::grid_file(std::string filename_, bool debugmode_)
:   filename(filename_), debugmode(debugmode_), datavalid(false),
    no_of_nodes(0), no_of_elements(0), nope(0), neta(0), nbou(0), nvel(0),
    min_x(-(std::numeric_limits<double>::max)()), 
    min_y(-(std::numeric_limits<double>::max)()), 
    min_depth(-(std::numeric_limits<double>::max)()), 
    max_x((std::numeric_limits<double>::max)()), 
    max_y((std::numeric_limits<double>::max)()), 
    max_depth((std::numeric_limits<double>::max)())
{
}

double grid_file::get_scale_factor(void)
{
    double i = 1.0;
    boost::uint32_t value = 0;
    double temp = (abs(max_x) > abs(max_y)) ? abs(max_x) : abs(max_y);

    value = boost::uint32_t(temp);

    while(value > 10)           // reducing everything to 1.something range
    {
        value = value / 10;
        i = i * 10;
    }
    return i;
}

///////////////////////////////////////////////////////////////////////////////
namespace util 
{
    template <typename T1>
    inline bool 
    read(std::istream &fort_14_file, T1& data1)
    {
        std::string instring;
        std::getline(fort_14_file, instring);
        if (fort_14_file.fail())
            return false;
            
        std::stringstream strstrm(instring);
        strstrm >> data1;
        return true;
    }

    template <typename T1, typename T2>
    inline bool 
    read(std::istream &fort_14_file, T1& data1, T2& data2)
    {
        std::string instring;
        std::getline(fort_14_file, instring);
        if (fort_14_file.fail())
            return false;
            
        std::stringstream strstrm(instring);
        strstrm >> data1 >> data2;
        return true;
    }

    template <typename T1, typename T2, typename T3>
    inline bool 
    read(std::istream &fort_14_file, T1& data1, T2& data2, T3& data3)
    {
        std::string instring;
        std::getline(fort_14_file, instring);
        if (fort_14_file.fail())
            return false;
            
        std::stringstream strstrm(instring);
        strstrm >> data1 >> data2 >> data3;
        return true;
    }

    template <typename T1, typename T2, typename T3, typename T4>
    inline bool 
    read(std::istream &fort_14_file, T1& data1, T2& data2, T3& data3, T4& data4)
    {
        std::string instring;
        std::getline(fort_14_file, instring);
        if (fort_14_file.fail())
            return false;
            
        std::stringstream strstrm(instring);
        strstrm >> data1 >> data2 >> data3 >> data4;
        return true;
    }

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    inline bool 
    read(std::istream &fort_14_file, T1& data1, T2& data2, T3& data3, T4& data4, 
        T5& data5)
    {
        std::string instring;
        std::getline(fort_14_file, instring);
        if (fort_14_file.fail())
            return false;
            
        std::stringstream strstrm(instring);
        strstrm >> data1 >> data2 >> data3 >> data4 >> data5;
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    void calc_minmax(std::vector<double>& vec, double& minvalue, double& maxvalue)
    {
        double minval = DBL_MAX;
        double maxval = -DBL_MAX;

        for (std::size_t i = 0; i < vec.size(); ++i) 
        {
            if (-99999 == vec[i])
                continue;
            
            minval = (std::min)(minval, vec[i]);
            maxval = (std::max)(maxval, vec[i]);
        }
    
        minvalue = minval;
        maxvalue = maxval;
    }
}

///////////////////////////////////////////////////////////////////////////////
// read data from input stream, read fort.14 file and populate data.
void grid_file::read(std::istream &fort_14_file, std::ostream* logfile)
{
    // Read in the grid name, no. of elements and no. of nodes
    if (!util::read(fort_14_file, grid_id, date) ||
        !util::read(fort_14_file, no_of_elements, no_of_nodes))
    {
        throw adcirc::exception("Unable to read from grid file: " + filename);
    }
    
    std::cerr << "Reading grid file: " << grid_id << " " << date << std::endl;
    if (logfile) 
        (*logfile) << "Reading grid file: " << grid_id << " " << date << std::endl;

    if (debugmode) {
        std::cerr << "elements: " << no_of_elements << std::endl;
        std::cerr << "nodes:    " << no_of_nodes << std::endl;
    }

    // loop to read in various x,y and depth value for all "no.of nodes"
    for (boost::uint32_t i = 0; i < no_of_nodes; ++i)
    {
        double temp, x_temp, y_temp, depth_temp;
        if (!util::read(fort_14_file, temp, x_temp, y_temp, depth_temp))
            throw adcirc::exception("Unable to read from grid file: " + filename);

        //  first element assign it to various variables (causes bugs when 
        //  negative values are present
        //  (thats why initialization to zero doesn't work)
        if (i == 0)
        {
            min_x = max_x = x_temp;
            min_y = max_y = y_temp;
            min_depth = max_depth = depth_temp;
        }

         x.push_back(x_temp);
         y.push_back(y_temp);
         depth.push_back(depth_temp);

         // calculate minimum and maximum values
         if (min_x > x_temp)
             min_x = x_temp;
         if (max_x < x_temp)
             max_x = x_temp;

         if (min_y > y_temp)
             min_y = y_temp;
         if (max_y < y_temp)
             max_y = y_temp;

         if (min_depth > depth_temp)
             min_depth = depth_temp;
         if (max_depth < depth_temp)
             max_depth = depth_temp;
    }

    if (debugmode) {
        std::cerr << "min vals: " << min_x << ", " << min_y << ", " << min_depth << std::endl;
        std::cerr << "max vals: " << max_x << ", " << max_y << ", " << max_depth << std::endl;
    }

    // read the connectivity information TOPOLOGY INFORMATION
    for (boost::uint32_t i = 0; i < no_of_elements; ++i)
    {
        boost::uint32_t temp = 0, nhy = 0, nm1 = 0, nm2 = 0, nm3 = 0;
        if (!util::read(fort_14_file, temp, nhy, nm1, nm2, nm3))
            throw adcirc::exception("Unable to read from grid file: " + filename);

//         NHY.push_back(nhy);
        NM1.push_back(nm1);
        NM2.push_back(nm2);
        NM3.push_back(nm3);
    }

    if (debugmode) {
        std::cerr << "topology information" << std::endl;
    }

    // boundary information
    
    // number of elevation specified boundary forcing segments.
    if (!util::read(fort_14_file, nope))
        throw adcirc::exception("Unable to read from grid file: " + filename);
    
    // total number of elevation specified boundary nodes
    if (!util::read(fort_14_file, neta))
        throw adcirc::exception("Unable to read from grid file: " + filename);

    if (debugmode) {
        std::cerr << "nope:     " << nope << std::endl;
        std::cerr << "neta:     " << neta << std::endl;
    }

    for (boost::uint32_t k = 0; k < nope; ++k) {
        boost::uint32_t nvdll_in = 0, ibtypee_in = 0;
        if (!util::read(fort_14_file, nvdll_in, ibtypee_in))
            throw adcirc::exception("Unable to read from grid file: " + filename);

        nvdll.push_back(nvdll_in);
        ibtypee.push_back(ibtypee_in);

        nbdv.push_back(std::vector<boost::uint32_t>());
        for (boost::uint32_t j = 0; j < nvdll_in; ++j) {
            boost::uint32_t nbdv_in = 0;
            if (!util::read(fort_14_file, nbdv_in))
                throw adcirc::exception("Unable to read from grid file: " + filename);
                
            nbdv.back().push_back(nbdv_in);
        }

        if (debugmode) {
            std::cerr << k << ": nvdll: " << nvdll_in << ", ibtypee: " 
                      << ibtypee_in << std::endl;
        }
    }
    
    // number of normal flow (discharge) specified boundary segments
    if (!util::read(fort_14_file, nbou))
        throw adcirc::exception("Unable to read from grid file: " + filename);
    
    nbvv.resize(nbou);
    barlanht.resize(nbou);
    barlancfsp.resize(nbou);
    ibconn.resize(nbou);
    barinht.resize(nbou);
    barincfsb.resize(nbou);
    barincfsp.resize(nbou);

    // total number of normal flow specified boundary nodes 
    if (!util::read(fort_14_file, nvel))
        throw adcirc::exception("Unable to read from grid file: " + filename);

    if (debugmode) {
        std::cerr << "nbou:     " << nbou << std::endl;
        std::cerr << "nvel:     " << nvel << std::endl;
    }

    for (boost::uint32_t k1 = 0; k1 < nbou; ++k1) {
        boost::uint32_t nvell_in = 0, ibtype_in = 0;
        if (!util::read(fort_14_file, nvell_in, ibtype_in))
            throw adcirc::exception("Unable to read from grid file: " + filename);
    
        nvell.push_back(nvell_in);
        ibtype.push_back(ibtype_in);
        
        nbvv[k1].resize(nvell_in);
        barlanht[k1].resize(nvell_in);
        barlancfsp[k1].resize(nvell_in);
        ibconn[k1].resize(nvell_in);
        barinht[k1].resize(nvell_in);
        barincfsb[k1].resize(nvell_in);
        barincfsp[k1].resize(nvell_in);

        for (boost::uint32_t j = 0; j < nvell_in; ++j) {
            boost::uint32_t nbvv_in = 0;
            double barlanht_in = 0, barlancfsp_in = 0;
            boost::uint32_t ibconn_in = 0;
            double barinht_in = 0, barincfsb_in = 0, barincfsp_in = 0;

            switch (ibtype_in) {
            case 0:  case 1:  case 2:
            case 10: case 11: case 12:
            case 20: case 21: case 22:
            case 30: case 31: case 32:
            case 40: case 41: case 42:
            case 50: case 51: case 52:
            case 60: case 61: case 62:
            case 70: case 71: case 72:
                if (!util::read(fort_14_file, nbvv_in))
                    throw adcirc::exception("Unable to read from grid file: " + filename);
                break;
                
            case 3:  case 13: case 23: 
            case 33: case 43: case 53: 
            case 63: case 73:
                if (!util::read(fort_14_file, nbvv_in, barlanht_in, barlancfsp_in))
                    throw adcirc::exception("Unable to read from grid file: " + filename);
                break;
                
            case 4:  case 14: case 24: 
            case 34: case 44: case 54: 
            case 64: case 74:
                if (!util::read(fort_14_file, nbvv_in, ibconn_in, barinht_in, barincfsb_in, barincfsp_in))
                    throw adcirc::exception("Unable to read from grid file: " + filename);
                break;
                
            default:
                throw adcirc::exception("Unknown ibtype: " + 
                    boost::lexical_cast<std::string>(ibtype_in));
            }
                
            nbvv[k1][j] = nbvv_in;
            barlanht[k1][j] = barlanht_in;
            barlancfsp[k1][j] = barlancfsp_in;
            ibconn[k1][j] = ibconn_in;
            barinht[k1][j] = barinht_in;
            barincfsb[k1][j] = barincfsb_in;
            barincfsp[k1][j] = barincfsp_in;
        }

        if (debugmode) {
            std::cerr << k1 << ": nvell: " << nvell_in << ", ibtype: " 
                      << ibtype_in << std::endl;
        }
    }
    datavalid = true;
}

///////////////////////////////////////////////////////////////////////////////
void grid_file::read (std::ostream* logfile)
{
    // handle plain or compressed ascii 
    std::string ext(boost::filesystem::extension(filename));
    if (ext != ".nc" && ext != ".cdf") {
        boost::iostreams::file_source fort_in(filename);
        boost::iostreams::filtering_istream fort_file;

        if (ext == ".gz") {
            fort_in.close();
            fort_in.open(filename, std::ios_base::in | std::ios_base::binary);
            fort_file.push(boost::iostreams::gzip_decompressor());
        }
        if(!fort_in.is_open())
            throw adcirc::exception("Unable to open the file: " + filename);

        fort_file.push(fort_in);
        read(fort_file, logfile);
        return;
    }

    // handle netcdf
     /*try {
        // open netCDF objects
       netCDF::NcFile ncfile(filename, netCDF::NcFile::read);

        // read dimensions
        no_of_nodes = (boost::uint32_t)read_file_dimension(ncfile, "node");
        no_of_elements = (boost::uint32_t)read_file_dimension(ncfile, "nele");

        // read attributes
        if (!read_file_attribute(ncfile, "grid", grid_id, false))
            read_file_attribute(ncfile, "model_domain", grid_id);
        date.clear();

        // read coordinates
        x.resize(no_of_nodes);
        y.resize(no_of_nodes);
        depth.resize(no_of_nodes);
        if (!read_netcdf_var(ncfile, "x", 1, 0, 0, 0, x, false))
            read_netcdf_var(ncfile, "lon", 1, 0, 0, 0, x);
        if (!read_netcdf_var(ncfile, "y", 1, 0, 0, 0, y, false))
            read_netcdf_var(ncfile, "lat", 1, 0, 0, 0, y);
        read_netcdf_var(ncfile, "depth", 1, 0, 0, 0, depth);

        // read grid topology
        NM1.resize(no_of_elements);
        NM2.resize(no_of_elements);
        NM3.resize(no_of_elements);
        if (read_netcdf_var(ncfile, "element", 2, 1, 0, 0, NM1, false)) {
            read_netcdf_var(ncfile, "element", 2, 1, 1, 0, NM2);
            read_netcdf_var(ncfile, "element", 2, 1, 2, 0, NM3);
        }
        else {
            read_netcdf_var(ncfile, "ele", 2, 1, 0, 0, NM1);
            read_netcdf_var(ncfile, "ele", 2, 1, 1, 0, NM2);
            read_netcdf_var(ncfile, "ele", 2, 1, 2, 0, NM3);
        }
        
        // calculate boundaries
        util::calc_minmax(x, min_x, max_x);
        util::calc_minmax(y, min_y, max_y);
        util::calc_minmax(depth, min_depth, max_depth);

        // set the rest to zero
        nope = 0;
        neta = 0;
        nbou = 0;
        nvel = 0;
    }
    catch (netCDF::exceptions::NcException const& e) {
        std::string msg("Caught NetCDF exception: ");
        msg += e.what();
        msg += " (file: " + filename + ")";
        if (logfile)
            (*logfile) << msg << std::endl;
        throw adcirc::exception(msg);
        }*/
}
///////////////////////////////////////////////////////////////////////////////
namespace util
{
    template <typename T1>
    bool write (std::ostream& outfile, T1 const& data1)
    {
        outfile << data1 << std::endl;
        return outfile.good();
    }

    template <typename T1, typename T2>
    bool write (std::ostream& outfile, T1 const& data1, T2 const& data2)
    {
        outfile << data1 << " " << data2 << std::endl;
        return outfile.good();
    }

    template <typename T1, typename T2, typename T3>
    bool write (std::ostream& outfile, T1 const& data1, T2 const& data2,
        T3 const& data3)
    {
        outfile << data1 << " " << data2 << " " << data3 << std::endl;
        return outfile.good();
    }

    template <typename T1, typename T2, typename T3, typename T4>
    bool write (std::ostream& outfile, T1 const& data1, T2 const& data2,
        T3 const& data3, T4 const& data4)
    {
        outfile << data1 << " " << data2 << " " << data3 << " " << data4 
                << std::endl;
        return outfile.good();
    }

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    bool write (std::ostream& outfile, T1 const& data1, T2 const& data2,
        T3 const& data3, T4 const& data4, T5 const& data5)
    {
        outfile << data1 << " " << data2 << " " << data3 << " " << data4 
                << " " << data5 << std::endl;
        return outfile.good();
    }

    template <typename T1, typename T2, typename T3, typename T4, typename T5, 
        typename T6>
    bool write (std::ostream& outfile, T1 const& data1, T2 const& data2,
        T3 const& data3, T4 const& data4, T5 const& data5, T6 const& data6)
    {
        outfile << data1 << " " << data2 << " " << data3 << " " << data4 
                << " " << data5 << " " << data6 << std::endl;
        return outfile.good();
    }
}

//  write the data to a file
void grid_file::write (std::ostream& outfile)
{
    if (!datavalid)
        throw adcirc::exception("Internal data is not valid!");
        
    if (!util::write (outfile, grid_id, date) ||
        !util::write (outfile, no_of_elements, no_of_nodes))
    {
        throw adcirc::exception("Unable to write to grid file: " + filename);
    }
    
    // loop to write out various x,y and depth value for all "no.of nodes"
    for (boost::uint32_t i = 0; i < no_of_nodes; ++i)
    {
        if (!util::write(outfile, i + 1, std::setprecision(8), x[i], y[i], depth[i]))
            throw adcirc::exception("Unable to write to grid file: " + filename);
    }

    // write the connectivity information 
    for (boost::uint32_t i = 0; i < no_of_elements; ++i)
    {
        if (!util::write(outfile, i + 1, /*NHY[i]*/3, NM1[i], NM2[i], NM3[i]))
            throw adcirc::exception("Unable to write to grid file: " + filename);
    }

    // number of elevation specified boundary forcing segments.
    if (!util::write(outfile, nope))
        throw adcirc::exception("Unable to write to grid file: " + filename);
    
    // total number of elevation specified boundary nodes
    if (!util::write(outfile, neta))
        throw adcirc::exception("Unable to write to grid file: " + filename);

    for (boost::uint32_t k = 0; k < nope; ++k) {
        if (!util::write(outfile, nvdll[k], ibtypee[k]))
            throw adcirc::exception("Unable to write to grid file: " + filename);
    
        for (boost::uint32_t j = 0; j < nvdll[k]; ++j) {
            if (!util::write(outfile, nbdv[k][j]))
                throw adcirc::exception("Unable to write to grid file: " + filename);
        }
    }
    
    // number of normal flow (discharge) specified boundary segments
    if (!util::write(outfile, nbou))
        throw adcirc::exception("Unable to write to grid file: " + filename);
    
    // total number of normal flow specified boundary nodes 
    if (!util::write(outfile, nvel))
        throw adcirc::exception("Unable to write to grid file: " + filename);

    for (boost::uint32_t k1 = 0; k1 < nbou; ++k1) {
        if (!util::write(outfile, nvell[k1], ibtype[k1]))
            throw adcirc::exception("Unable to write to grid file: " + filename);
    
        for (boost::uint32_t j = 0; j < nvell[k1]; ++j) {
            switch (ibtype[k1]) {
            case 0:  case 1:  case 2:
            case 10: case 11: case 12:
            case 20: case 21: case 22:
            case 30: case 31: case 32:
            case 40: case 41: case 42:
            case 50: case 51: case 52:
            case 60: case 61: case 62:
            case 70: case 71: case 72:
                if (!util::write(outfile, nbvv[k1][j]))
                    throw adcirc::exception("Unable to write to grid file: " + filename);
                break;
                
            case 3:  case 13: case 23: 
            case 33: case 43: case 53: 
            case 63: case 73:
                if (!util::write(outfile, nbvv[k1][j], std::setprecision(8), 
                    barlanht[k1][j], barlancfsp[k1][j]))
                {
                    throw adcirc::exception("Unable to write to grid file: " + filename);
                }
                break;
                
            case 4:  case 14: case 24: 
            case 34: case 44: case 54: 
            case 64: case 74:
                if (!util::write(outfile, nbvv[k1][j], ibconn[k1][j], 
                    std::setprecision(8), barinht[k1][j], barincfsb[k1][j], 
                    barincfsp[k1][j]))
                {
                    throw adcirc::exception("Unable to write to grid file: " + filename);
                }
                break;
                
            default:
                throw adcirc::exception("Unknown ibtype: " + 
                    boost::lexical_cast<std::string>(ibtype[k1]));
            }
        }
    }
}

void grid_file::write (std::ostream* logfile)
{
    // handle plain or compressed ascii 
    std::string ext(boost::filesystem::extension(filename));
    if (ext != ".nc" && ext != ".cdf") {
        boost::iostreams::file_sink fort_out(filename);
        boost::iostreams::filtering_ostream fort_file;

        if (ext == ".gz") {
            fort_out.close();
            fort_out.open(filename, std::ios_base::out | std::ios_base::binary);
            fort_file.push(boost::iostreams::gzip_compressor());
        }
        if(!fort_out.is_open())
            throw adcirc::exception("Unable to open the file: " + filename);

        fort_file.push(fort_out);
        write(fort_file);
        return;
    }

    // write new NetCDF file
    /*try {
        // create netCDF objects
        netCDF::NcFile ncfile(filename, 
            netCDF::NcFile::FileMode(netCDF::NcFile::write|netCDF::NcFile::replace),
            netCDF::NcFile::classic64);

//         // read dimensions
//         no_of_nodes = util::read_file_dimension(ncfile, "node");
//         no_of_elements = util::read_file_dimension(ncfile, "nele");
// 
//         // read attributes
//         util::read_file_attribute(ncfile, "grid", grid_id);
//         date.clear();
// 
//         // read coordinates
//         util::read_netcdf_var(ncfile, "x", 1, 0, x);
//         util::read_netcdf_var(ncfile, "y", 1, 0, y);
//         util::read_netcdf_var(ncfile, "depth", 1, 0, depth);
// 
//         // read grid topology
//         util::read_netcdf_var(ncfile, "element", 2, 0, NM1);
//         util::read_netcdf_var(ncfile, "element", 2, 1, NM2);
//         util::read_netcdf_var(ncfile, "element", 2, 2, NM3);
    }
    catch (netCDF::exceptions::NcException const& e) {
        std::string msg("Caught NetCDF exception: ");
        msg += e.what();
        msg += " (file: " + filename + ")";
        if (logfile)
          (*logfile) << msg << std::endl;
        throw adcirc::exception(msg);
    }*/
}

///////////////////////////////////////////////////////////////////////////////
}   // namespace adcirc

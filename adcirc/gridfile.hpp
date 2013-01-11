//     Copyright (c) 2006 Hartmut Kaiser. Distributed under the Boost
//     Software License, Version 1.0. (See accompanying file
//     LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef GRID_FILE_HPP_INCLUDED
#define GRID_FILE_HPP_INCLUDED

#include <vector>
#include <string>
#include <iosfwd>

#include <boost/cstdint.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace adcirc 
{
    ///////////////////////////////////////////////////////////////////////////
    // Class to read the grid file
    struct grid_file
    {
        std::string filename;
        bool debugmode;
        bool datavalid;
                
        // Various data read from fort.14 file
        std::string grid_id;          // Grid name
        std::string date;             // grid file version
        boost::uint32_t no_of_nodes;   
        boost::uint32_t no_of_elements;

        // Geometry information
        // std::vector index corresponds to the node number
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> depth;

        // Topology information
        // std::vector index corresponds to the junction element (JE) number 
        // NHY = number of nodes per element. 
        //    At present the only allowable value for the number of nodes per 
        //    element is 3 indicating a triangular element with linear basis 
        //    functions.
//         std::vector<boost::uint32_t> NHY;

        //  NM_1, NM_2, NM_3 = node numbers comprising element JE. 
        //  These must be specified in a counter clockwise direction moving 
        //  around the element.
        std::vector<int> NM1;
        std::vector<int> NM2;
        std::vector<int> NM3;

        // boundary nodes
        unsigned int nope;                            // number of elevation specified boundary forcing segments.
        unsigned int neta;                            // total number of elevation specified boundary nodes

        std::vector<boost::uint32_t> nvdll;              // number of nodes in elevation boundary segment k
        std::vector<boost::uint32_t> ibtypee;            // elevation boundary type
        std::vector<std::vector<boost::uint32_t> > nbdv; // node numbers on elevation specified boundary segment k.

        // normal flow boundary segments
        unsigned int nbou;                            // number of normal flow (discharge) specified boundary segments
        unsigned int nvel;                            // total number of normal flow specified boundary nodes 

        std::vector<boost::uint32_t> nvell;              // number of nodes in normal flow specified boundary 
        std::vector<boost::uint32_t> ibtype;             // boundary type

        std::vector<std::vector<boost::uint32_t> > nbvv; // node numbers on normal flow boundary segment
        std::vector<std::vector<double> > barlanht;   // external barrier height
        std::vector<std::vector<double> > barlancfsp; // coefficient of free surface supercritical flow at external barrier node
        std::vector<std::vector<boost::uint32_t> > ibconn;  // back face node paired with the front face node
        std::vector<std::vector<double> > barinht;    // internal barrier height 
        std::vector<std::vector<double> > barincfsb;  // coefficient of free surface sub-critical flow at internal barrier node
        std::vector<std::vector<double> > barincfsp;  // coefficient of free surface supercritical flow at internal barrier node 

        // Minimum and Maximum values for various vectors
        double min_x, min_y, min_depth, max_x, max_y, max_depth;

        double get_scale_factor(void);

        // Read function to read fort.14 file and populate data
        grid_file(std::string filename, bool debugmode_ = false);

        void read (std::istream& infile, std::ostream* logfile = NULL);
        void write (std::ostream& outfile);

        // better read/write allowing to use netcdf 
        void read (std::ostream* logfile = NULL);
        void write (std::ostream* logfile = NULL);
    }; 

///////////////////////////////////////////////////////////////////////////////
}   // end of namespace adcirc

#endif  // !GRID_FILE_HPP_INCLUDED

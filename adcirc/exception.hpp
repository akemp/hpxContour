//     Copyright (c) 2006-2009 Hartmut Kaiser. Distributed under the Boost
//     Software License, Version 1.0. (See accompanying file
//     LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef EXCEPTION_HPP_INCLUDED
#define EXCEPTION_HPP_INCLUDED

#include <exception>
#include <string>

///////////////////////////////////////////////////////////////////////////////
namespace adcirc 
{
    ///////////////////////////////////////////////////////////////////////////
    class exception : public std::exception
    {
    private:
        std::string msg;
        
    public:
        exception (std::string m) : msg(m) 
        {
        }
        ~exception() throw() {}
        
        virtual char const *what() const throw() { return msg.c_str(); }
    };

///////////////////////////////////////////////////////////////////////////////
}   // namespace adcirc

#endif

#ifndef HEADER
#define HEADER

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "adcirc/gridfile.hpp"
#include "adcirc/outputfile.hpp"

#include <boost/cstdint.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/extensions/index/rtree/rtree.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

using namespace boost::geometry;
using namespace std;

typedef model::point<double, 2, cs::cartesian> pointxy;
typedef std::pair<pointxy, pointxy> line;
typedef model::box<pointxy> box;

static bool point2Equals(pointxy p1, pointxy p2)
{
    return ((get<0>(p1) == get<0>(p2)) && (get<1>(p1) == get<1>(p2)));
}

struct triangle
{
    triangle(vector<double> input)
    {
        if (input.size() != 9)
        {
            model::point<double, 3, cs::cartesian> p1;
            assign_values(p1, 0.0, 0.0, 0.0);
            model::point<double, 3, cs::cartesian> p2;
            assign_values(p2, 0.0, 1.0, 0.0);
            model::point<double, 3, cs::cartesian> p3;
            assign_values(p3, 1.0, 0.0, 0.0);
            model::point<double, 3, cs::cartesian> p4;
            assign_values(p4, 0.0, 0.0, 0.0);
            this->points.push_back(p1);
            this->points.push_back(p2);
            this->points.push_back(p3);
            this->points.push_back(p4);
            append(this->poly,this->points);
            cout << "Error: Invalid input. Using default value.\n";
        }
        else
        {
            model::point<double, 3, cs::cartesian> p1;
            assign_values(p1, input[0], input[1], input[2]);
            model::point<double, 3, cs::cartesian> p2;
            assign_values(p2, input[3], input[4], input[5]);
            model::point<double, 3, cs::cartesian> p3;
            assign_values(p3, input[6], input[7], input[8]);
            model::point<double, 3, cs::cartesian> p4;
            assign_values(p4, input[0], input[1], input[2]);
            this->points.push_back(p1);
            this->points.push_back(p2);
            this->points.push_back(p3);
            this->points.push_back(p4);
            append(this->poly,this->points);
        }
        Zaverage = 0;
        maximum = -10000;
        minimum = 10000;
        for (int i = 0; i < 3; i++)
        {
            float temp = get<2>(points[i]);
            Zaverage += temp;
            if (temp > maximum)
                maximum = temp;
            if (temp < minimum)
                minimum = temp;
            
        }
        Zaverage /= 3.0;
        
        
    }
    bool intersect(double hei)
    {
        if (minimum <= hei && maximum >= hei)
        {
            return true;
        }
        return false;
    }
    pointxy intRatio (model::point<double, 3, cs::cartesian> p1, model::point<double, 3, cs::cartesian> p2, double height, int& intersects)
    {
        double z1 = get<2>(p1), z2 = get<2>(p2);
        intersects = 1;
        pointxy seg;
        if (z1 >= height && z2 <= height)
        {
            if (z1 == z2)
                assign_values(seg, get<0>(p1), get<1>(p1));
            else
            {
                double range = z1-z2;
                double diff1 = z1-height;
                double ratio1 = 1-diff1/range;
                double diff2 = height-z2;
                double ratio2 = 1-diff2/range;
                double x = get<0>(p1)*ratio1+get<0>(p2)*(ratio2);
                double y = get<1>(p1)*ratio1+get<1>(p2)*(ratio2);
                assign_values(seg, x, y);
            }
        }
        else if (z2 >= height && z1 <= height)
        {
            if (z1 == z2)
                assign_values(seg, get<0>(p1), get<1>(p1));
            else
            {
                double range = z2-z1;
                double diff1 = z2-height;
                double ratio1 = 1-diff1/range;
                double diff2 = height-z1;
                double ratio2 = 1-diff2/range;
                double x = get<0>(p1)*ratio2+get<0>(p2)*(ratio1);
                double y = get<1>(p1)*ratio2+get<1>(p2)*(ratio1);
                assign_values(seg, x, y);
            }
        }
        else
        {
            assign_values(seg, 0, 0);
            intersects = 0;
        }
        return seg;
    }
    vector<pointxy> planeInter(double height)
    {
        model::point<double, 3, cs::cartesian> p1 = points[0];
        model::point<double, 3, cs::cartesian> p2 = points[1];
        model::point<double, 3, cs::cartesian> p3 = points[2];

        vector<pointxy > s;
        int intersects = 0;
        pointxy seg = intRatio(p1,p2,height,intersects);
        if (intersects)
            s.push_back(seg);
        seg = intRatio(p2,p3,height,intersects);
        if (intersects)
            s.push_back(seg);
        seg = intRatio(p3,p1,height,intersects);
        if (intersects)
            s.push_back(seg);
        vector<pointxy > final;
        for (int i = 0; i < s.size(); i++)
        {
            bool contains = false;
            pointxy cur = s[i];
            for (int j = i+1; j < s.size(); j++)
            {
                if (point2Equals(s[j],cur))
                {
                    contains = true;
                    break;
                }
            }
            if (!contains)
                final.push_back(cur);
        }
        return final;
    }

    static bool compare(const triangle& i, const triangle& j)
    {
        return (i.minimum < j.minimum);
    }

    model::polygon<model::point<double, 3, cs::cartesian> > poly;
    vector<model::point<double, 3, cs::cartesian>> points;
    double Zaverage;
    double maximum;
    double minimum;
};

struct sliver
{
    sliver(double height, double minX, double maxX, double minY, double maxY, int place){
        this->height = height;
        this->place = place;
        this->indexer = 0;
    }

    vector<pointxy> intersect(triangle trig, double height)
    {
        vector<pointxy> ints = trig.planeInter(height);
        if (ints.size() == 2)
        {
            intersects.push_back(std::make_pair(ints[0],ints[1]));
        }
        else if (ints.size() == 3)
        {
            intersects.push_back(std::make_pair(ints[0],ints[1]));
            intersects.push_back(std::make_pair(ints[1],ints[2]));
            intersects.push_back(std::make_pair(ints[2],ints[0]));
        }
        indexer++;
        return ints;
    }

    double height;
    //vector<vector<pointxy>> intersects;
    vector<line> intersects;
    int place;
    int indexer;
};
#endif

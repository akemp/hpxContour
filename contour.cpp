#include <math.h>

#ifndef INCLUDED
#define INCLUDED
#include "datastructures.hpp"
#endif

#include <boost/accumulators/accumulators.hpp>

#include <boost/geometry/extensions/index/rtree/rtree.hpp>

//using namespace spatial_index;

typedef model::point<double, 2, cs::cartesian> pointxy;
typedef vector<pointxy> lineS;
typedef std::pair<pointxy, pointxy> line;
typedef model::box<pointxy> box;
typedef index::rtree<box, int> rtree;


void runPseudoRandomTests(int testNum);
void createTestData(int x, int y);
void runDataTests(double iterations, int limit);

//Find the line segments which share vertices and dynamically create an R-Tree to sort through them.
//Remove and add vertices as appropriate
void findLine(vector<pointxy>& ints, const double &margin, rtree& rt, int& counts, vector<lineS> &segs)
{
    pointxy p1(ints[0].get<0>()-margin,ints[0].get<1>()-margin);
    pointxy p2(ints[0].get<0>()+margin,ints[0].get<1>()+margin);
    box b1(p1,p2);
    pointxy p3(ints[1].get<0>()-margin,ints[1].get<1>()-margin);
    pointxy p4(ints[1].get<0>()+margin,ints[1].get<1>()+margin);
    box b2(p3,p4);
    auto v = rt.find(b1);
    //if **NO** shared vertex with the first point found
    if (v.size() <= 0)
    {
        v = rt.find(b2);
        //if **NO** shared vertex with the second point found
        if (v.size() <= 0)
        {
            if (boost::geometry::distance(ints[0], ints[1]) > margin)
            {
                rt.insert(b1, counts);
                rt.insert(b2, counts);
                lineS temp;
                temp.push_back(ints[0]);
                temp.push_back(ints[1]);
                segs.push_back(temp);
                counts++;
            }
        }
        else
        {
            int index = v[0];
            //remove the shared point, we have a new endpoint now
            rt.remove(b2,index);
            //add the new endpoint
            rt.insert(b1,index);
            pointxy begin = segs[index].front();
            pointxy end = segs[index].back();
            double dist1 = boost::geometry::distance(begin, ints[1]);
            double dist2 = boost::geometry::distance(end, ints[1]);
            if (dist2 < dist1)
            {
                //A little error checking. Due to rounding errors this might be toggled a sporadically during an application.
                double dist = boost::geometry::distance(end, ints[1]);
                if (dist > margin)
                    cout << "Error 2-1: Distance too great!\n";
                segs[index].push_back(ints[0]);
            }
            else
            {
                reverse(segs[index].begin(), segs[index].end());
                end = segs[index].back();
                double dist = boost::geometry::distance(end, ints[1]);
                //A little error checking. Due to rounding errors this might be toggled a sporadically during an application.
                if (dist > margin)
                {
                    cout << "Error 2-2: Distance too great!\n";
                }
                segs[index].push_back(ints[0]);
            }
        }
    }
    else
    {
        int index = v[0];
        //remove the shared point, we have a new endpoint now
        rt.remove(b1,index);
        //add the new endpoint
        rt.insert(b2,index);
        pointxy begin = segs[index].front();
        pointxy end = segs[index].back();
        double dist1 = boost::geometry::distance(begin, ints[0]);
        double dist2 = boost::geometry::distance(end, ints[0]);
        if (dist2 < dist1)
        {
            double dist = boost::geometry::distance(end, ints[0]);
            //A little error checking. Due to rounding errors this might be toggled a sporadically during an application.
            if (dist > margin)
                cout << "Error 1-1: Distance too great!\n";
            segs[index].push_back(ints[1]);
        }
        else
        {
            reverse(segs[index].begin(), segs[index].end());
            end = segs[index].back();
            double dist = boost::geometry::distance(end, ints[0]);
            //A little error checking. Due to rounding errors this might be toggled a sporadically during an application.
            if (dist > margin)
                cout << "Error 1-2: Distance too great!\n";
            segs[index].push_back(ints[1]);
        }
    }
}

template <typename num_t>
void outputLines(const size_t numLines)
{
    // read the triangular grid fro the given file
    adcirc::grid_file gridfile("ec95d_grid.14");
    gridfile.read();

    // now you can access the coordinates for the triangles
    boost::uint32_t num_nodes = gridfile.no_of_nodes;         // number of nodes
    std::vector<double> x_coords = gridfile.x;                // X coordinates
    std::vector<double> y_coords = gridfile.y;                // Y coordinates
    std::vector<double> z_coords = gridfile.depth;                // Y coordinates

    // these 3 vectors hold the node numbers for the triangles
    boost::uint32_t num_triangles = gridfile.no_of_elements;  // number of triangles 
    std::vector<int> nodes1 = gridfile.NM1;
    std::vector<int> nodes2 = gridfile.NM2;
    std::vector<int> nodes3 = gridfile.NM3;
        
    double maxX = -100000;
    double minX = 100000;
    double maxY = -100000;
    double minY = 100000;
    double maxZ = -100000;
    double minZ = 100000;
    vector<triangle> trigs;
    //We fill up the geometry containers and find the range of data here. We need the range for when we generate the contour lines.
    for (boost::uint32_t i = 0; i < num_triangles; ++i)
    {
        vector<double> dummy;
        int node1 = nodes1[i]-1;
        int node2 = nodes2[i]-1;
        int node3 = nodes3[i]-1;
        double x = x_coords[node1];
        double y = y_coords[node1];
        double z = z_coords[node1];
        dummy.push_back(x);
        dummy.push_back(y);
        dummy.push_back(z);
        if (z < minZ)
            minZ = z;
        if (z > maxZ)
            maxZ = z;
        x = x_coords[node2];
        y = y_coords[node2];
        z = z_coords[node2];
        dummy.push_back(x);
        dummy.push_back(y);
        dummy.push_back(z);
        if (z < minZ)
            minZ = z;
        if (z > maxZ)
            maxZ = z;
        x = x_coords[node3];
        y = y_coords[node3];
        z = z_coords[node3];
        dummy.push_back(x);
        dummy.push_back(y);
        dummy.push_back(z);
        if (z < minZ)
            minZ = z;
        if (z > maxZ)
            maxZ = z;
        triangle trig(dummy);
        trigs.push_back(trig);
    }//*/
    
    double inc = (maxZ-minZ)/numLines;

    double inv = 1.0/inc;
    
    vector<sliver> slivers;

    double place = minZ;

    vector<rtree> rts;
    vector<int> counts;
    vector<vector<lineS>> manySegs;
    

    //fill up the various storage containers with empty members
    for (double i = 0; i < numLines; ++i)
    {
        rtree rt(32,8);
        rts.push_back(rt);
        slivers.push_back(sliver(place,minX,maxX,minY,maxY,i) );
        place += inc;
        counts.push_back(0);
        vector<lineS> temp;
        manySegs.push_back(temp);
    }
    
    double margin = 0.000001; //margin of error
    
    //Go through the triangles and generate the contour lines.
    //This could be parallelized by running ALL of the triangles through several different threads which only generate contour lines for
    //a set height. For example, 4 threads with each one processing a certain quartile of data.
    for (vector<triangle>::iterator it = trigs.begin(); it < trigs.end(); ++it)
    {
        triangle trig = (*it);
        for (int i = min((int)(floor((trig.minimum - minZ)*inv)),0); i <= (int)(ceil((trig.maximum - minZ)*inv)) && i < numLines; ++i)
        {
            double height = slivers[i].height;
            if (trig.minimum <= slivers[i].height && trig.maximum >= slivers[i].height)
            {
                vector<pointxy> ints = slivers[i].intersect(trig,height);
                if (ints.size() == 2)
                {
                    findLine(ints,margin,rts[i],counts[i],manySegs[i]);
                }
                else if (ints.size() == 3)
                {
                    //Planar intersections - just create three separate line segments and deal with those.
                    cout << "Planar intersection located.\n";
                    vector<pointxy> ints1;
                    vector<pointxy> ints2;
                    vector<pointxy> ints3;
                    ints1.push_back(ints[0]);
                    ints1.push_back(ints[1]);

                    ints2.push_back(ints[1]);
                    ints2.push_back(ints[2]);

                    ints3.push_back(ints[2]);
                    ints3.push_back(ints[0]);

                    findLine(ints1,margin,rts[i],counts[i],manySegs[i]);
                    findLine(ints2,margin,rts[i],counts[i],manySegs[i]);
                    findLine(ints3,margin,rts[i],counts[i],manySegs[i]);
                }
            }
        }
    }
    
    cout << "Complete.\n";

    for (int i = 0; i < manySegs.size(); i++)
    {
        cout << "Entry " << i+1 << endl;
        bool found = true;
        bool begun = true;
        //This is almost a carbon copy of findLine. The only difference is that it uses a conditional so it does not stop
        //combining line segments with shared vertices prematurely. Accuracy is key here.
        while (found || begun)
        {
            int counter = 0;
            begun = false;
            found = false;
            rtree rt(32,8);
            vector<lineS> segs;
            for (int j = 0; j < manySegs[i].size(); j++)
            {
                lineS temp = manySegs[i][j];
                pointxy p1(temp.front().get<0>()-margin,temp.front().get<1>()-margin);
                pointxy p2(temp.front().get<0>()+margin,temp.front().get<1>()+margin);
                box b1(p1,p2);
                pointxy p3(temp.back().get<0>()-margin,temp.back().get<1>()-margin);
                pointxy p4(temp.back().get<0>()+margin,temp.back().get<1>()+margin);
                box b2(p3,p4);
                auto v = rt.find(b1);
                if (v.size() <= 0)
                {
                    v = rt.find(b2);
                    if (v.size() <= 0)
                    {
                        rt.insert(b1, counter);
                        rt.insert(b2, counter);
                        segs.push_back(temp);
                        counter++;
                    }
                    else
                    {
                        found = true;
                        int index = v[0];
                        rt.remove(b2,index);
                        rt.insert(b1,index);
                        lineS holder = segs[index];
                        pointxy begin = holder.front();
                        pointxy end = holder.back();
                        double dist1 = boost::geometry::distance(begin, temp.back());
                        double dist2 = boost::geometry::distance(end, temp.back());
                        if (dist2 < dist1)
                        {
                            reverse(temp.begin(), temp.end());
                            double dist = boost::geometry::distance(holder.back(), temp.front());
                            if (dist > margin)
                                cout << "Error 2-1: Distance too great!\n";
                            holder.insert(holder.end(), temp.begin()+1, temp.end());

                        }
                        else
                        {
                            reverse(temp.begin(), temp.end());
                            reverse(holder.begin(), holder.end());
                            double dist = boost::geometry::distance(holder.back(), temp.front());
                            if (dist > margin)
                            {
                                cout << "Error 2-2: Distance too great!\n";
                            }
                            holder.insert(holder.end(), temp.begin()+1, temp.end());
                        }
                        segs[index] = holder;
                    }
                }
                else
                {
                    found = true;
                    int index = v[0];
                    rt.remove(b1,index);
                    rt.insert(b2,index);
                    lineS holder = segs[index];
                    pointxy begin = holder.front();
                    pointxy end = holder.back();
                    double dist1 = boost::geometry::distance(begin, temp.front());
                    double dist2 = boost::geometry::distance(end, temp.front());
                    if (dist2 < dist1)
                    {
                        double dist = boost::geometry::distance(holder.back(), temp.front());
                        if (dist > margin)
                            cout << "Error 1-1: Distance too great!\n";
                        holder.insert(holder.end(), temp.begin()+1, temp.end());
                    }
                    else
                    {
                        reverse(holder.begin(), holder.end());
                        double dist = boost::geometry::distance(holder.back(), temp.front());
                        if (dist > margin)
                            cout << "Error 1-2: Distance too great!\n";
                        holder.insert(holder.end(), temp.begin()+1, temp.end());
                    }
                    segs[index] = holder;
                }
            }
            if (found)
                manySegs[i] = segs;
        }
    }
    //we output the file to a simple graphic format. The segment vectors are not grouped by height in manySegs
    ofstream fout("out.obj");

    int count = 1;

    for (int i = 0; i < slivers.size(); i++)
    {
        double height = -slivers[i].height/200.0;
        int j = 0;
        for (vector<lineS>::iterator it = manySegs[i].begin(); it < manySegs[i].end(); ++it)
        {
            lineS intersect = (*it);
            for (int k = 0; k < intersect.size(); k++)
            {
                fout << "v " << intersect[k].get<0>() << " " << intersect[k].get<1>() << " " << height << endl;
            }
            fout << "g g_" << count << endl;
            for (int k = 0; k < intersect.size()-1; k++)
            {
                fout << "l " << count << " ";
                count++;
                fout << count << endl;
            }
            count++;
        }
    }
    fout.close();
    
}

int main()
{
    try {
        
        outputLines<double>(110);
        return 0;
    }
    catch (adcirc::exception const& e) {
        std::cerr << "contour: caught exception: " << e.what() << std::endl;
    }
}

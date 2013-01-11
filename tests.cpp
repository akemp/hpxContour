#ifndef INCLUDED
#define INCLUDED

#include "datastructures.hpp"

#endif


bool compare(const triangle& i, const triangle& j);
void createTestData(int x, int y)
{
    ofstream fout("test1.txt");
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            fout << i << " ";
        }
        fout << endl;
    }
    fout.close();

    fout.open("test2.txt");
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            fout << j << " ";
        }
        fout << endl;
    }
    fout.close();

    fout.open("test3.txt");
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            fout << (i+j)%2 << " ";
        }
        fout << endl;
    }
    fout.close();
}
vector<model::point<double, 2, cs::cartesian>> intersect(triangle trig, int height)
{
    vector<model::point<double, 2, cs::cartesian>> s = trig.planeInter(height);
    return s;
}//*/

void runPseudoRandomTests(int testNum)
{
    bool err = false;
    cout << "Nonplanar intersection - 1 side higher than others.\n";
    for (int i = 1; i < testNum; i++)
    {
        int height = rand()%100000;
        vector<double> dummy;
        dummy.push_back(i);//1X
        dummy.push_back(i);//1Y
        dummy.push_back(height+i*2);//1Z

        dummy.push_back(i);//2X
        dummy.push_back(i+1);//2Y
        dummy.push_back(height+i*2+1);//2Z

        dummy.push_back(i+1);//3X
        dummy.push_back(i+1);//3Y
        dummy.push_back(height+i*2);//3Z
        triangle trig(dummy);

        vector<model::point<double, 2, cs::cartesian>> intersects = intersect(trig,height+i*2);
        if (intersects.size() != 2)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 2!\n";
            err = true;
            break;
        }
        intersects = intersect(trig,height+i*2+1);
        if (intersects.size() != 1)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 1!\n";
            err = true;
            break;
        }
    }
        
    cout << "Planer intersection - no side higher than others.\n";
    for (int i = 1; i < testNum; i++)
    {
        int height = rand()%100000;
        vector<double> dummy;
        dummy.push_back(i);//1X
        dummy.push_back(i);//1Y
        dummy.push_back(i*2+1+height);//1Z

        dummy.push_back(i);//2X
        dummy.push_back(i+1);//2Y
        dummy.push_back(i*2+1+height);//2Z

        dummy.push_back(i+1);//3X
        dummy.push_back(i+1);//3Y
        dummy.push_back(i*2+1+height);//3Z
        triangle trig(dummy);

        vector<model::point<double, 2, cs::cartesian>> intersects = intersect(trig,i*2+height);
        if (intersects.size() != 0)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 0!\n";
            err = true;
            break;
        }
        intersects = intersect(trig,i*2+1+height);
        if (intersects.size() != 3)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 3!\n";
            err = true;
            break;
        }
    }
        
    cout << "Skewed intersection - sides have different heights.\n";
    for (int i = 1; i < testNum; i++)
    {
        int height = rand()%100000;
        vector<double> dummy;
        dummy.push_back(i);//1X
        dummy.push_back(i);//1Y
        dummy.push_back(height+i*2-1);//1Z

        dummy.push_back(i);//2X
        dummy.push_back(i+1);//2Y
        dummy.push_back(height+i*2);//2Z

        dummy.push_back(i+1);//3X
        dummy.push_back(i+1);//3Y
        dummy.push_back(height+i*2+1);//3Z
        triangle trig(dummy);

        vector<model::point<double, 2, cs::cartesian>> intersects = intersect(trig,i*2-1+height);
        if (intersects.size() != 1)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 1!\n";
            err = true;
            break;
        }
        intersects = intersect(trig,i*2+height);
        if (intersects.size() != 2)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 2!\n";
            err = true;
            break;
        }
        intersects = intersect(trig,i*2+1+height);
        if (intersects.size() != 1)
        {
            cout << "Iteration: " << i << endl;
            cout << "Error! Intersects total " << intersects.size() << " when they should be 1!\n";
            err = true;
            break;
        }
    }
    if (!err)
        cout << "All pseudorandom tests completed without errors.\n";
}

struct vertex {
    vertex(int x,int y, int z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    } 
    bool vertex::operator==(const vertex &rhs) {
        
        return (rhs.x == this->x && rhs.y == this->y && rhs.z == this->z); 
  }
    int x, y, z;
};
void runDataTests(double iterations, int limit)
{
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
        
    double maxZ = -100000;
    double minZ = 100000;
    vector<vector<triangle>> trigArray;
    vector<vertex> combinations;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vertex temp1((i+j)%3,(i+j+1)%3,(i+j+2)%3);
            bool contains = false;
            for (int k = 0; k < combinations.size(); k++)
            {
                vertex temp = combinations[k];
                if (temp1 == temp)
                {
                    contains = true;
                    break;
                }
            }
            if (!contains)
            {
                combinations.push_back(temp1);
            }
            vertex temp2((i+j+2)%3,(i+j)%3,(i+j+1)%3);
            contains = false;
            for (int k = 0; k < combinations.size(); k++)
            {
                vertex temp = combinations[k];
                if (temp == temp2)
                {
                    contains = true;
                    break;
                }
            }
            if (!contains)
            {
                combinations.push_back(temp2);
            }
            vertex temp3((i+j+2)%3,(i+j+1)%3,(i+j)%3);
            contains = false;
            for (int k = 0; k < combinations.size(); k++)
            {
                vertex temp = combinations[k];
                if (temp == temp3)
                {
                    contains = true;
                    break;
                }
            }
            if (!contains)
            {
                combinations.push_back(temp3);
            }
        }
    }
    vector<vector<triangle>> trigAr;
    {
        for (int j = 0; j < combinations.size(); j++)
        {
            vector<triangle> trigs;
            // this will print the coordinates for all triangles
            for (boost::uint32_t i = 0; i < num_triangles && (i < limit || limit <= 0); ++i)
            {
                vector<double> dummy;
                vector<int> nodes;
                nodes.push_back(nodes1[i]-1);
                nodes.push_back(nodes2[i]-1);
                nodes.push_back(nodes3[i]-1);
                double z = z_coords[nodes[combinations[j].x]];
                dummy.push_back(x_coords[nodes[combinations[j].x]]);
                dummy.push_back(y_coords[nodes[combinations[j].x]]);
                dummy.push_back(z);
                if (z > maxZ)
                    maxZ = z;
                if (z < minZ)
                    minZ = z;
                z = z_coords[nodes[combinations[j].y]];
                dummy.push_back(x_coords[nodes[combinations[j].y]]);
                dummy.push_back(y_coords[nodes[combinations[j].y]]);
                dummy.push_back(z);
                z = z_coords[nodes[combinations[j].z]];
                if (z > maxZ)
                    maxZ = z;
                if (z < minZ)
                    minZ = z;
                dummy.push_back(x_coords[nodes[combinations[j].z]]);
                dummy.push_back(y_coords[nodes[combinations[j].z]]);
                dummy.push_back(z);
                if (z > maxZ)
                    maxZ = z;
                if (z < minZ)
                    minZ = z;
                triangle trig(dummy);
                trigs.push_back(trig);

            }//*/
            trigAr.push_back(trigs);
        }
    }
    for (int i = 0; i < trigAr.size(); i++)
        sort(trigAr[i].begin(), trigAr[i].end(), triangle::compare);
    int count = 0;
    for (int z = 0; z < trigAr.size(); z++)
    {
        for (vector<triangle>::iterator it = trigAr[z].begin(); it < trigAr[z].end(); ++it)
        {
            vector<model::point<double, 3, cs::cartesian> > points = (*it).points;
            for (int j = 0; j < 3; j++)
            {
                if (get<2>(points[j]) < minZ)
                    minZ = get<2>(points[j]);
                if (get<2>(points[j]) > maxZ)
                    maxZ = get<2>(points[j]);
            }
        }
    }
    bool fail = false;
    vector<vector<vector<model::point<double, 2, cs::cartesian>>>> matches;
    for (int z = 0; z < trigAr.size(); z++)
    {
        {
            vector<vector<model::point<double, 2, cs::cartesian>>> temp;
            matches.push_back(temp);
        }
        for (double height = minZ; height < maxZ; height += (maxZ-minZ)/iterations)
        {
            for (vector<triangle>::iterator it = trigAr[z].begin(); it < trigAr[z].end(); ++it)
            {
                triangle trig = (*it);
                vector<model::point<double, 2, cs::cartesian>> intersects = intersect(trig,height);
                for (int i = 0; i < intersects.size(); ++i)
                {
                    model::point<double, 2, cs::cartesian> cur = intersects[i];
                    for (int j = i+1; j < intersects.size(); ++j)
                    {
                        if (equals(cur, intersects[j]))
                        {
                            fail = true;
                            break;
                        }
                    }
                }
                matches[z].push_back(intersects);
            }
        }
    }
    
    if (fail)
    {
        cout << "Error! Equal points recorded!\n";
        return;
    }
    else
    {
        cout << "Tests completed without incident.\n";
    }
    for (int i = 0; i < matches.size(); i++)
    {
        for (int j = 0; j < matches.size(); j++)
        {
            if (matches[i].size() != matches[j].size())
            {
                cout << "Error! Number of matches for permuted version " << i << " not equal to permuted version " << j << "!\n";
                return;
            }
        }
    }
    for (int j = 0; j < matches.size(); j++)
    {
        for (int k = 0; k < matches.size(); ++k)
        {
            for (int z = 0; z < matches[j].size(); ++z)
            {
                if ((matches[k][z]).size() != (matches[j][z]).size())
                {
                    cout << (matches[j][z]).size() << endl;
                    cout << (matches[k][z]).size() << endl;
                    cout << "Error! Number of matches for permuted version " << j << " not equal to permuted version " << k << "!\n";
                }
            }
        }
    }
}

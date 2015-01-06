#include "datastructures.hpp"
#include <shapefil.h>
#include <boost/algorithm/string.hpp>
#include <netcdf.h>

bool objFlag()
{
	return true;
}
bool shpFlag()
{
	return true;
}
bool svgFlag()
{
	return false;
}


bool inputValues(variables_map vm, string &input, string &output, vector<string> &attnames, vector<double> &heights,
	int &timestep, bool &outputInput, string& xs, string& ys)
{
	/*
  options_description cmdline("usage: " HPX_APPLICATION_STRING " [options]");
  
  cmdline.add_options()
    ( "input"
      , value<string>()
      , "Input filename")

    ( "output"
      , value<string>()->default_value("output")
      , "Output filename prefix")

    ( "atts"
      , value<string>()->default_value("depth")
      , "Attributes to read in")
          
      ( "heights"
      , value<string>()->default_value("50,60,70,80")
      , "Heights to generate the contours from")

      ( "timestep"
      , value<int>()->default_value(0)
      , "Timestep to run from")
    ;
	*/
	
	timestep = 0;
	output = "out";
	if (!vm.count("input"))
	{
		cout << "Error: No input specified!\n";
		return false;
	}
	input = vm["input"].as<string>();
	if (vm.count("timestep"))
		timestep = vm["timestep"].as<int>();
	if (vm.count("output"))
		output = vm["output"].as<string>();
	string atts = "depth";
	if (vm.count("atts"))
		atts = vm["atts"].as<string>();
	xs = "x";
	if (vm.count("x"))
		xs = vm["x"].as<string>();
	ys = "y";
	if (vm.count("y"))
		ys = vm["y"].as<string>();
	vector<string> holder;
	boost::split(holder, atts, boost::is_any_of(","), boost::token_compress_on);
	for (int i = 0; i < holder.size(); ++i)
	{
		attnames.push_back(holder[i]);
	}
	if (attnames.size() <= 0)
	{
		cout << "Invalid attribute inputs.\n";
		return false;
	}
	string height ="50,100,150,200,250,300,350,400,450,500";
	if (vm.count("heights"))
		height = vm["heights"].as<string>();
	
	if (height.find(".txt") != string::npos)
	{
		std::ifstream fin(height);
		while (fin.good())
		{
			double temp = 0;
			fin >> temp;
			heights.push_back(temp);
		}
		fin.close();
	}
	else
	{
		holder.clear();
		boost::split(holder, height, boost::is_any_of(","), boost::token_compress_on);
		for (int i = 0; i < holder.size(); ++i)
		{
			double temp = atof(holder[i].c_str());
			if (temp != 0)
				heights.push_back(temp);
		}
	}
	if (vm["outputInput"].as<string>().length() > 1)
	{
		outputInput = true;
	}
    std::sort(heights.begin(), heights.end());
	return true;
}

bool readNC(boost::shared_array<double> &x, boost::shared_array<double> &y, vector<boost::shared_array<double>> &depth, boost::shared_array<int> &ele,
	int &node, int &nele, string input, vector<string> attnames, int timestep, bool outputInput, string xs, string ys)
{
   /* This will be the netCDF ID for the file and data variable. */
     int status, ncid, ndims, nvars, ngatts, unlimdimid;

   /* Loop indexes, and error handling. */
   int retval;
   int time;
   bool hasTime = false;
   if ((retval = nc_open(input.c_str(), NC_NOWRITE, &ncid)))
	   return false;
     nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
	 for (int i = 0; i < ndims; ++i)
	 {
		 using namespace std;
		 char* name = new char[1024];
		 size_t lenp = 0;
		 nc_inq_dim(ncid,i,name,&lenp);
		 string finder = name;
		 if (finder.find("nele") != -1)
			 nele = lenp;
		 else if (finder.find("node") != -1)
			 node = lenp;
		 else if (finder.find("time") != -1)
		 {
			 time = lenp;
			 if (time > 1)
			 hasTime = true;
		 }
	 }
	 

   /* Get the varid of the data variable, based on its name. */

	 x = boost::shared_array<double>(new double[node]);
	 y = boost::shared_array<double>(new double[node]);
	 for (int i = 0; i < attnames.size(); ++i)
		 depth.push_back(boost::shared_array<double>(new double[node]));
	 ele = boost::shared_array<int>(new int[nele*3]);
	 
	 for (int i = 0; i < nvars; ++i)
	 {
		 using namespace std;
		 char* name = new char[1024];
		 nc_inq_varname(ncid,i,name);
		 string finder = name;
		 bool found = false;
		 {
			for (int j = 0; j < attnames.size(); ++j)
			{ 
				if (finder.find(attnames[j]) != -1)
				{
					if (attnames[j] == "depth" || !hasTime)
					{	size_t end = node;
						size_t start = 0;
						nc_get_vara_double (ncid, i, &start, &end, depth[j].get());
						found = true;
					}
					else
					{
						size_t end[] = {1,node};
						size_t start[] = {timestep,0};
						nc_get_vara_double (ncid, i, start, end, depth[j].get());
						found = true;
					}

				}
			}
		 }
		 if (!found)
		 {
			 if (finder.find(xs) != -1 && finder.size() == xs.size())
			 {
				 size_t end = node;
				 size_t start = 0;
				 nc_get_vara_double (ncid, i, &start, &end, x.get());
			 }
			 else if (finder.find(ys) != -1 && finder.size() == xs.size())
			 {
				 size_t end = node;
				 size_t start = 0;
				 nc_get_vara_double (ncid, i, &start, &end, y.get());
			 }
			 else if (finder.find("element") != -1)
			 {
				 size_t end[] = {nele,3};
				 size_t start[] = {0,0};
				 nc_get_vara_int (ncid, i, start, end, ele.get());
			 }
		 }
	 }
   if ((retval = nc_close(ncid)))
	   return false;
   if (outputInput)
   {   
	string filename = "debOutput.obj";
	std::ofstream fout(filename);
	int count = 1;
	cout << "Writing to file.\n";
	int counter = 0;
	for (int i = 1; i < nele/3; i++)
	{
		int index = i*3;
		int ind = ele.get()[index]-1;
		double x1 = x.get()[ind];
		double y1 = y.get()[ind];
		double z1 = depth[0].get()[ind];
		index++;
		ind = ele.get()[index]-1;
		double x2 = x.get()[ind];
		double y2 = y.get()[ind];
		double z2 = depth[0].get()[ind];
		index++;
		ind = ele.get()[index]-1;
		double x3 = x.get()[ind];
		double y3 = y.get()[ind];
		double z3 = depth[0].get()[ind];
		if (z1 > -999 && z2 > -999 && z3 > -999)
		{
			fout << "v " << x1 << " " << z1 << " " << y1 << endl;
			fout << "v " << x2 << " " << z2 << " " << y2 << endl;
			fout << "v " << x3 << " " << z3 << " " << y3 << endl;
			fout << "f " << count << " " << count + 1 << " " << count + 2 << endl;
			count += 3;
		}
	}
	fout.close();
	cout << "Written to file.\n";
   }
   return true;
}

void outputObj(boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> &contours, vector<double> &heights, int nL,string output)
{
	string filename = output + ".obj";
	std::ofstream fout(filename);
	int count = 1;
	cout << "Writing to file.\n";
	int counter = 0;
	for (int lin = 0; lin < nL; ++lin)
	{
		for (vector<vector<pointxy>>::iterator it = contours[lin]->begin(); it < contours[lin]->end(); ++it)
		{
			vector<pointxy> intersect = (*it);
			for (int k = 0; k < intersect.size(); k++)
			{
				fout << "v " << intersect[k].get<0>() << " " << -heights[counter] << " " << intersect[k].get<1>() << endl;
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
		counter++;
	}
	fout.close();
	cout << "Written to file.\n";
}
void outputShape(boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> &contours, vector<double> &heights, int nL, string output)
{
	vector<vector<double>> inds;
	{
		SHPHandle	hSHPHandle;
		string filename = output + ".shp";
		hSHPHandle = SHPCreate( filename.c_str(), SHPT_POLYGON );
	
		int counter = 0;
		int total = 0;
		for (int lin = 0; lin < nL; ++lin)
		{
			vector<double> ind;
			int count = 0;
			for (vector<vector<pointxy>>::iterator it = contours[lin]->begin(); it < contours[lin]->end(); ++it)
			{
				SHPObject	*psShape;
				double	xa[10000], ya[10000], za[10000], ma[10000];
				vector<pointxy> intersect = (*it);
				int size = 0;
				for (int k = 0; k < intersect.size(); ++k)
				{
					xa[k] = intersect[k].get<0>();
					ya[k] = intersect[k].get<1>();
					za[k] = 0;
					ma[k] = 0;
					size++;
				}
				psShape = SHPCreateObject( SHPT_POLYGON, -1, 0, NULL, NULL,
											size, xa, ya, za, ma );
				SHPWriteObject( hSHPHandle, -1, psShape );
				SHPDestroyObject( psShape );
				count++;
				ind.push_back(total);
				total++;
			}
			counter++;
			inds.push_back(ind);
		}
		SHPClose( hSHPHandle );
	}
	{
		DBFHandle	hDBF;
		string filename = output + ".dbf";
		hDBF = DBFCreate(filename.c_str());
		DBFAddField( hDBF, "height", FTDouble, 10, 10 );
		for (int i = 0; i < heights.size(); ++i)
		{
			double height = heights[i];
			for (int j = 0; j < inds[i].size(); ++j)
			{
				DBFWriteDoubleAttribute(hDBF, inds[i][j], 0, height);
			}
		}
		DBFClose(hDBF);
	}
	return;
}

void outputSvg(boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> &contours,
    vector<double> &heights, int nL, string output)
{
     {
        vector<string> colors;
        colors.push_back("black");
        colors.push_back("gray");
        colors.push_back("dimgray");
        colors.push_back("lightgray");
	    string filename = output + ".svg";
	    std::ofstream fout(filename);
	    cout << "Writing to file.\n";
        double multer = 100;
        double divx = 1.0;
        double divy = 1.0;
        fout << " <svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"10000\" height = \"10000\">\n";
	    int counter = 0;
        //fout << "<rect width=\"" << multer << "\" height=\"" << multer << "\" style=\"fill:Sienna;fill-opacity=0.5;stroke-width:0;\"/>";
        double addy = 0;
        double addx = 10000;
        {
	        for (int lin = 0; lin < nL; ++lin)
	        {
		        for (vector<vector<pointxy>>::iterator it = contours[lin]->begin(); it < contours[lin]->end(); ++it)
		        {
                    int random = rand()%colors.size();
                    int random2 = rand()%colors.size();
			        vector<pointxy> intersect = (*it);
                    {
			            fout << "<path d=\"";
                        fout << "M" << intersect[0].get<0>()*multer/divx+addx << "," << intersect[0].get<1>()*multer/divy+addy << " L";
			            for (int k = 0; k < intersect.size(); k++)
			            {
				            fout << " " << intersect[k].get<0>()*multer/divx+addx << "," << intersect[k].get<1>()*multer/divy+addy << " ";
			            }
                        fout << " " << intersect.back().get<0>()*multer/divx+addx << "," << intersect.back().get<1>()*multer/divy+addy << " ";
                        fout << " \" stroke=\"";
                        fout << colors[random];
                        fout << "\" stroke-width=\"2\" stroke-opacity=\"" << 1 << "\" fill=\"" << "none";

                        fout << "\" fill-opacity=\"" << 0 << "\"/>" << endl;
                    }
                    counter++;
		        }
	        }
        }
        fout << "</svg>";
	    fout.close();
	    cout << "Written to svg file.\n";
     }
}

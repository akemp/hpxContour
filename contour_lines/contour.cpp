#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/lcos.hpp>

#include "datastructures.hpp"

bool objFlag();
bool shpFlag();
bool svgFlag();
bool inputValues(variables_map vm, string &input, string &output, vector<string> &attnames, vector<double> &heights,
	int &timestep, bool &outputInput, string& xs, string& ys);

void outputObj(boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> &contours, vector<double> &heights, int nL,string output);
void outputShape(boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> &contours, vector<double> &heights, int nL,string output);
void outputSvg(boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> &contours,
    vector<double> &heights, int nL, string output);

int retLine(double current, int size,
	boost::shared_array<double> &x, boost::shared_array<double> &y, vector<boost::shared_array<double>> &depth,
	boost::shared_array<int> &ele, boost::shared_ptr< vector <vector <pointxy> >>& contours);

bool readNC(boost::shared_array<double> &x, boost::shared_array<double> &y, vector<boost::shared_array<double>> &depth, boost::shared_array<int> &ele,
	int &node, int &nele, string input, vector<string> attnames, int timestep, bool outputInput, string xs, string ys);


void outputLines(string input, string output, vector<string> attnames, vector<double> heights, int timestep, bool outputInput, string xs, string ys)
{
    size_t const os_threads = hpx::get_os_thread_count();
    cout << "Program will run on " << os_threads << " threads.\n";
	int nL = heights.size();
	cout << "Reading in file.\n";


   /* Open the file. NC_NOWRITE tells netCDF we want read-only access
    * to the file.*/
	 
	boost::shared_array<double> x;
	boost::shared_array<double> y;
	vector<boost::shared_array<double>> depth;
	boost::shared_array<int> ele;
	 int node = 0;
	 int nele = 0;
	 if (!readNC(x,y,depth,ele,node,nele,input, attnames, timestep, outputInput,xs,ys))
	 {
		 cout << "File read error. Check input filename and validity of input data.\n";
		 return;
	 }
    cout << "File read. Starting timer and generating contour lines.\n";
	hpx::util::high_resolution_timer t;
		
	boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>> contours =
		boost::shared_array<boost::shared_ptr<vector<vector<pointxy>>>>(new boost::shared_ptr<vector<vector<pointxy>>>[nL]);
   {
		using hpx::lcos::future;
		using hpx::async;
		vector<future<int>> waiting;
	   for (int i = 0; i < heights.size(); ++i)
	   {
		   /*
int retLine(double current, int size, int counter,
	boost::shared_array<double> &x, boost::shared_array<double> &y, boost::shared_array<double> &depth,
	boost::shared_array<int> &ele, boost::shared_array<vector<vector<pointxy>>> &contours)
		   */
		   contours[i] = boost::shared_ptr<vector<vector<pointxy>>>(new vector<vector<pointxy>>());
		   waiting.push_back(async(&retLine,heights[i], nele/3, x,y,depth,ele,contours[i]));
	  }
	   hpx::wait_all(waiting);
   }
	double ti = t.elapsed();
	cout << "Program completed in " << ti << "s.\n";
	std::ofstream fout("times.txt", std::ios_base::app);
	fout << os_threads << " " << ti << endl;
	fout.close();
	if (objFlag())
		outputObj(contours, heights, nL,output);
	if (shpFlag())
		outputShape(contours, heights, nL,output);
	if (svgFlag())
		outputSvg(contours, heights, nL,output);
	return;
}

int hpx_main(variables_map& vm)
{
	string input, output;
	vector<string> attnames;
	vector<double> heights;
	string xs, ys;
	int timestep;
  // Configure application-specific options.
	bool outputInput = false;;
	if (!inputValues(vm, input, output, attnames, heights, timestep, outputInput, xs, ys))
	{
		hpx::finalize(); // Handles HPX shutdown
		return 1;
	}
	outputLines(input,output,attnames,heights,timestep, outputInput,xs, ys);
	hpx::finalize(); // Handles HPX shutdown
	return 0;
}
int main(int argc, char *argv[])
{
	options_description desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");
	desc_commandline.add_options()
    ("input", value<string>(), "Set input file\n")
    ("output", value<string>()->default_value("out"), "Set output file\n")
    ("x", value<string>()->default_value("x"), "Set x coord attribute name. Default: x\n")
    ("y", value<string>()->default_value("y"), "Set y coord attribute name. Default: y\n")
    ("atts", value<string>()->default_value("depth"), "Set attributes to read. Mulitple attributes are separated by a comma (,)\n")
    ("heights", value<string>()->default_value("50,100,150,200,250,300,350,400,450,500"), "Heights to generate the contours from. Mulitple heights are separated by a comma (,)\n")
    ("timestep", value<int>()->default_value(0), "Which timestep to run (if avaliable)\n")
    ("outputInput", value<string>()->default_value(""), "Whether to output the input file as an obj file - \"true\" for yes (don't use this if you don't want to output the input variables)\n")
;
	return hpx::init(desc_commandline, argc, argv);
}

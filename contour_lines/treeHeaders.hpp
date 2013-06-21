#ifndef HEADER
#define HEADER

#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

#include <boost/shared_array.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

using std::vector;
using std::cout;
using std::endl;
using std::string;

using namespace boost::geometry;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)
typedef model::point<double, 2, cs::cartesian> pointxy;
#endif
struct kd_tree
{
  struct node
  {
    int i0, i1, split;

    node() : i0(-1), i1(-1), split(-1)
    {}
    friend bool is_leaf(const node &n) { return n.split==int(-1); }
  };
  pointxy min, max;
  std::vector<node> nodes;

  kd_tree(const vector<pointxy> &points);
  int nearest(const pointxy &p,
           const vector<pointxy> &points,
           int index = -1,
		   double distance = 100000000) const;
};
struct build_task;

struct lower_x;

struct lower_y;
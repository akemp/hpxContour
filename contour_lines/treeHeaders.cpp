#ifndef HEADER
#define HEADER
#define NOMINMAX
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

using std::cout;
using std::endl;
using std::vector;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(boost::geometry::cs::cartesian)
typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> pointxy;

#endif

inline int next_pow_of_2(int n)
{
  --n;  //the following works for 32 bit ints
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return ++n;
}


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
  vector<node> nodes;

  kd_tree(const vector<pointxy> &points);
  int nearest(const pointxy &p,
           const vector<pointxy> &points,
           int index = -1,
           double distance = 100000000) const;
};

struct build_task
{
  int first, last, node, dim;
};

inline void expand(pointxy &min, pointxy &max, const pointxy &p)
{
  if (p.get<0>()<min.get<0>()) min.set<0>(p.get<0>());
  else if (p.get<0>()>max.get<0>()) max.set<0>(p.get<0>());
  if (p.get<1>()<min.get<1>()) min.set<1>(p.get<1>());
  else if (p.get<1>()>max.get<1>()) max.set<1>(p.get<1>());
  return;
}

struct lower_x
{
  const std::vector<pointxy> &points;

  lower_x(const std::vector<pointxy> &p) : points(p)
  {}
  bool operator()(int i1, int i2) const
  { return points[i1].get<0>()<points[i2].get<0>(); }
};

struct lower_y
{
  const std::vector<pointxy> &points;

  lower_y(const std::vector<pointxy> &p) : points(p)
  {}
  bool operator()(int i1, int i2) const
  { return points[i1].get<1>()<points[i2].get<1>(); }
};

#define maxS(a,b)            (((a) > (b)) ? (a) : (b))

#define minS(a,b)            (((a) < (b)) ? (a) : (b))


kd_tree::kd_tree(const std::vector<pointxy> &points)
{
  const int stack_size = 32;
  const int count = points.size();
  if (count==0)
    return;
  int *index = new int[count];
  for (int i = 0; i<count; ++i)
    index[i] = i;
  const int m = next_pow_of_2(count);
  const int node_count = minS( int(m - 1), int(2*count - m/2 - 1) );
  nodes.reserve(node_count);
  nodes.push_back(node());
  build_task tasks[stack_size] = {{/*first*/0, /*last*/count - 1, /*node*/0, /*dim*/0}};
  int current_task = 0;
  pointxy min(100000, 100000);
  pointxy max(-100000, -100000);
  do
  {
    build_task task = tasks[current_task];
    node &n = nodes[task.node];
    if (task.last - task.first<=1)
    {
      n.i0 = index[task.first];
      expand(min, max, points[n.i0]);
      if (task.first!=task.last)
      {
        n.i1 = index[task.last];
        expand(min, max, points[n.i0]);
      }
      assert(is_leaf(n));
      --current_task;
      continue;
    }
    const int k = (task.first + task.last)/2;
    if (task.dim==0)
      std::nth_element(index + task.first, index + k, index + task.last + 1, lower_x(points));
    else
      std::nth_element(index + task.first, index + k, index + task.last + 1, lower_y(points));
    const int
      i0 = nodes.size(),
      i1 = i0 + 1;
    n.split = index[k];
    n.i0 = i0;
    n.i1 = i1;
    const int next_dir = task.dim^1;
    const build_task
      task0 = {task.first, k, i0, next_dir},
      task1 = {k + 1, task.last, i1, next_dir};
    tasks[current_task] = task0;
    nodes.push_back(node());
    ++current_task;
    assert(current_task<stack_size);
    tasks[current_task] = task1;
    nodes.push_back(node());
  } while (current_task!=-1);
  delete [] index;
  assert(nodes.size()==node_count);
  this->min = min;
  this->max = max;
}

struct traverse_task
{
  pointxy c;
  double dist;
  int node, dim;
};

inline double bind_to_segment(double x, double x0, double x1)
{
  if (x<x0)
    return x0;
  else if (x>x1)
    return x1;
  return x;
}

inline pointxy bind_to_rect(const pointxy &p, const pointxy &min, const pointxy &max)
{
  return pointxy(bind_to_segment(p.get<0>(), min.get<0>(), max.get<0>()), bind_to_segment(p.get<1>(), min.get<1>(), max.get<1>()));
}

int kd_tree::nearest(const pointxy &p,
                        const vector<pointxy> &points,
                        int index,
                        double distance) const
{
  const int stack_size = 32;
  const pointxy c = bind_to_rect(p, min, max);
  traverse_task tasks[stack_size] = {{c, boost::geometry::distance(p, c), /*node*/0, /*dim*/0}};
  int current_task = 0;
  int best = -1;
  do
  {
    if (tasks[current_task].dist>distance)
      continue;
    for (;;)
    {
      const traverse_task &task = tasks[current_task];
      const node &n = nodes[task.node];
      if (is_leaf(n))
      {
        assert(n.i0!=int(-1));
		const double d = boost::geometry::distance(p, points[n.i0]);
        if (d<distance && n.i0 > index)
        {
          distance = d;
          best = n.i0;
        }
        if (n.i1!=int(-1) && n.i1 > index)
        {
          const double d = boost::geometry::distance(p, points[n.i1]);
          if (d<distance)
          {
            distance = d;
            best = n.i1;
          }
        }
        break;
      }
      traverse_task nearest, farther;
      pointxy c = task.c;
      bool lower_is_nearest;
      if (task.dim==0)
      {
        const double split = points[n.split].get<0>();
        c.set<0>(split);
        lower_is_nearest = p.get<0>()<=split;
      }
      else
      {
        const double split = points[n.split].get<1>();
        c.set<1>(split);
        lower_is_nearest = p.get<1>()<=split;
      }
      if (lower_is_nearest)
      {
        nearest.node = n.i0;
        farther.node = n.i1;
      }
      else
      {
        nearest.node = n.i1;
        farther.node = n.i0;
      }
      const int next_dim = task.dim^1;
      nearest.c = task.c;
      nearest.dist = 0.;
      nearest.dim = next_dim;
      farther.c = c;
	  farther.dist = boost::geometry::distance(p, farther.c);
      farther.dim = next_dim;
      tasks[current_task] = farther;
      ++current_task;
      assert(current_task<stack_size);
      tasks[current_task] = nearest;
    }
  } while (--current_task!=-1);
  return best;
}

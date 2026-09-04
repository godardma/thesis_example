#include <codac>
using namespace codac2;
using namespace std;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Boolean_set_operations_2.h> 
#include <CGAL/squared_distance_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::Point_2 Point;

typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;

typedef CGAL::Exact_predicates_tag Itag;

typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDTBase;
typedef CGAL::Constrained_triangulation_plus_2<CDTBase> CDT;

typedef CDT::Finite_faces_iterator Face_iterator;
typedef CDT::Vertex_handle Vertex_handle;

typedef CGAL::Polygon_set_2<Kernel> Polygon_set_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;


typedef CGAL::Partition_traits_2<Kernel>                         Traits;
typedef Traits::Polygon_2                                   Polygon_2_T;
typedef Traits::Point_2                                     Point_2;
typedef std::list<Polygon_2_T>                                Polygon_list;

typedef CGAL::Gps_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

Polygon to_polygon (const Parallelepiped& p)
{
    return Polygon({p.vertices()[0],p.vertices()[1],p.vertices()[3],p.vertices()[2]});
}

vector<Polygon> to_polygon (const vector<Parallelepiped>& v_par)
{
    vector<Polygon> polygons;
    for (auto& p : v_par)
        polygons.push_back(to_polygon(p));
    return polygons;
}

Polygon_2 oriented_polygon(const Polygon_2& polygon)
{
  Polygon_2 oriented_polygon = polygon;
  if (oriented_polygon.is_clockwise_oriented())
    oriented_polygon.reverse_orientation();
  return oriented_polygon;
}

Point to_point(const Vector& v) {
    return Point(v[0], v[1]);
}

Vector to_vector(const Point& p) {
    double p0 = CGAL::to_double(p.x());
    double p1 = CGAL::to_double(p.y());
    return Vector({p0,p1});
}

Parallelepiped to_parallelepiped (const IntervalVector& x)
{
    Vector z ({x[0].mid(), x[1].mid()});
    Matrix A({{(x[0].diam())/2.0, 0},
                  {0, x[1].diam()/2.0}});
    return Parallelepiped (z, A);
}

bool are_almost_adjacent(const Polygon_2& p1, const Polygon_2& p2, const double epsilon = 1e-6) {
    // 1. First, use do_intersect() to check for true adjacency or overlap
    for (auto e1 = p1.edges_begin(); e1 != p1.edges_end(); ++e1) {
        for (auto e2 = p2.edges_begin(); e2 != p2.edges_end(); ++e2) {
            if (CGAL::do_intersect(*e1, *e2)) {
                return true;
            }
        }
    }

    // 2. Compute the minimum distance between all edge pairs
    Kernel::FT min_sq_dist = std::numeric_limits<double>::max();  // or FT::max() if using FT
    for (auto e1 = p1.edges_begin(); e1 != p1.edges_end(); ++e1) {
        for (auto e2 = p2.edges_begin(); e2 != p2.edges_end(); ++e2) {
            Kernel::FT sq_dist = CGAL::squared_distance(*e1, *e2);
            if (sq_dist < min_sq_dist) {
                min_sq_dist = sq_dist;
            }
        }
    }

    // 3. Compare with epsilon
    return (CGAL::to_double(min_sq_dist) <= epsilon * epsilon);
}

Polygon_2 to_cgal(const Parallelepiped& p) 
{
    Polygon_2 polygon;

    auto vertices = p.vertices();
    polygon.push_back(to_point(vertices[0]));
    polygon.push_back(to_point(vertices[1]));
    polygon.push_back(to_point(vertices[3]));
    polygon.push_back(to_point(vertices[2]));

    return oriented_polygon(polygon);
}

Polygon_2 to_cgal (const IntervalVector& x)
{
    Polygon_2 polygon;

    polygon.push_back(Point(x[0].lb(), x[1].lb()));
    polygon.push_back(Point(x[0].ub(), x[1].lb()));
    polygon.push_back(Point(x[0].ub(), x[1].ub()));
    polygon.push_back(Point(x[0].lb(), x[1].ub()));


    return oriented_polygon(polygon);
}

Polygon_2 to_cgal (const Polygon& polygon)
{
    Polygon_2 cgal_polygon;

    for (const auto& p : polygon.vertices()) 
        cgal_polygon.push_back(Point(CGAL::to_double(p[0].mid()), CGAL::to_double(p[1].mid())));

    return oriented_polygon(cgal_polygon);
}

vector<Vector> to_codac (const Polygon_2& polygon)
{
  std::vector<Vector> points;

  for (const Point& p : polygon) 
    points.push_back(to_vector(p));

  return points;
}

Polygon to_codac (const Face_iterator& it)
{
  Point v0 = it->vertex(0)->point();
  Point v1 = it->vertex(1)->point();
  Point v2 = it->vertex(2)->point();

  Polygon polygon({
      to_vector(v0),
      to_vector(v1),
      to_vector(v2)
  });

  return polygon;
}

Polygon_list get_convex_partition (const Polygon_2& polygon)
{
  Polygon_list  partition_polys;
  CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                              polygon.vertices_end(),
                              std::back_inserter(partition_polys));
  return partition_polys;
}

vector<Point> get_points (const Polygon_2& polygon)
{
  vector<Point> points;
  for (auto vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit)
    points.push_back(*vit);
  return points;
}

Polygon_2 convert_polygon (const Polygon_2_T& polygon)
{
  Polygon_2 new_polygon;
  for (const Point& p : polygon) 
    new_polygon.push_back(Point(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
  return new_polygon;
}

vector<Polygon> get_triangle (const CGAL::Constrained_Delaunay_triangulation_2<Kernel>& cdt)
{
  vector<Polygon> v_polygons;
  for (Face_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
    v_polygons.push_back(to_codac(it));
  return v_polygons;
}

Polygon_with_holes_2 generate_polygon_with_hole (const vector<Parallelepiped>& v_par)
{
  Polygon_set_2 polygon_set;

  int i = 0;

  for (auto& p : v_par)
    polygon_set.join(to_cgal(p));
  

  std::list<Polygon_with_holes_2> result;
  polygon_set.polygons_with_holes(std::back_inserter(result));
  return *result.begin();
}

Polygon_with_holes_2 generate_polygon_with_hole (const vector<IntervalVector>& v_b)
{
  Polygon_set_2 polygon_set;

  int i = 0;

  for (auto& b : v_b)
    polygon_set.join(to_cgal(b));
  

  std::list<Polygon_with_holes_2> result;
  polygon_set.polygons_with_holes(std::back_inserter(result));
  return *result.begin();
}

double winding_angle (const Vector &p, const vector<Parallelepiped> v_par)
{
  double angle = 0;
  for (int i = 0; i < v_par.size(); ++i)
  {
    Vector z1 = v_par[i].c;
    Vector z2 = v_par[(i+1)%v_par.size()].c;

    double d_ang = (atan2(z2[1]-p[1], z2[0]-p[0])-atan2(z1[1]-p[1], z1[0]-p[0]));
    if (d_ang > M_PI)
      d_ang -= 2*M_PI;
    else if (d_ang < -M_PI)
      d_ang += 2*M_PI;
    angle += d_ang;
  }
  return angle;
}

bool is_inside (const vector<Polygon>& polygons, const vector<Parallelepiped> v_par)
{
  auto p =(((polygons[0]).vertices())[0]).mid();
  int winding_number = std::round(winding_angle(p, v_par) / (2*M_PI));
  return winding_number!=0;
}

vector<vector<Polygon>> triangulate_holes(const vector<Parallelepiped>& v_par)
{
  vector<vector<Polygon>> hole_connex_polygons;


  Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_par);

  int i = 0;
  
  for (auto it = poly_with_holes.holes_begin(); it != poly_with_holes.holes_end(); ++it) 
  {
    Polygon_2 hole = oriented_polygon(*it);
    
    Polygon_list  partition_polys = get_convex_partition(hole);

    vector<Polygon> polygons;
    for (auto it2 = partition_polys.begin(); it2 != partition_polys.end(); ++it2)
    {
        auto convex_polygon = *it2;
        CGAL::Constrained_Delaunay_triangulation_2<Kernel> cdt;

        std::vector<Vector> vertices=to_codac(convert_polygon(convex_polygon));
        std::vector<Point> points=get_points(convert_polygon(convex_polygon));

        cdt.insert(points.begin(), points.end());

        for (auto& p : get_triangle (cdt))
            polygons.push_back(p);
    }
    hole_connex_polygons.push_back(polygons);
  }

  return hole_connex_polygons;
}

bool is_on_bounded_side_or_boundary (const Point& p, const Polygon_2& polygon)
{

    if (polygon.bounded_side(p) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(p) == CGAL::ON_BOUNDARY)
        return true;
    return false;
}

vector<Polygon> triangulate_outside (const vector<Parallelepiped>& v_par, const IntervalVector& init_box)
{
  vector<Polygon> exterior_polygons;

  Polygon_2 init_poly = to_cgal(init_box);

  Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_par);

  Polygon_2 outer_boundary = oriented_polygon(poly_with_holes.outer_boundary());

  CDT cdt_box;
  for (std::size_t i = 0; i < init_poly.size(); ++i)
      cdt_box.insert_constraint(init_poly[i], init_poly[(i+1)%init_poly.size()]);
  
  for (std::size_t i = 0; i < outer_boundary.size(); ++i)
      cdt_box.insert_constraint(outer_boundary[i], outer_boundary[(i+1)%outer_boundary.size()]);


  // Polygon_list  partition_polys_outer = get_convex_partition(outer_boundary);

  for (auto fit = cdt_box.finite_faces_begin(); fit != cdt_box.finite_faces_end(); ++fit) 
  {

    Point a = fit->vertex(0)->point();
    Point b = fit->vertex(1)->point();
    Point c = fit->vertex(2)->point();
    // m is the center of the triangle
    Point m = CGAL::centroid(a, b, c);

    bool point_a_is_inside = is_on_bounded_side_or_boundary(a, outer_boundary);
    bool point_b_is_inside = is_on_bounded_side_or_boundary(b, outer_boundary);
    bool point_c_is_inside = is_on_bounded_side_or_boundary(c, outer_boundary);
    bool point_m_is_inside = is_on_bounded_side_or_boundary(m, outer_boundary);


    // if no points are inside a polygon, draw the triangle
    if (!point_a_is_inside || !point_b_is_inside || !point_c_is_inside || !point_m_is_inside)
        exterior_polygons.push_back(to_codac(fit));
    }

  return exterior_polygons;
}
    
vector<Parallelepiped> eliminate_fake_parallelepipeds(const vector<Parallelepiped>& v_par, const vector<vector<Polygon>>& outside_polygons, vector<vector<Polygon>>& inside_polygons)
{
  vector<Parallelepiped> v_par_no_fake;
  vector<Parallelepiped> v_par_fake;

  for (auto& p : v_par)
  {
    Polygon_2 poly_p = to_cgal(p);
    bool is_fake = true;
    for (const auto& polygon : outside_polygons)
    {
      for (const auto& triangle : polygon)
      {
        Polygon_2 poly_q = to_cgal(triangle);
        if (are_almost_adjacent(poly_p, poly_q))
        {
          is_fake = false;
          break;
        }
      }
      if (!is_fake)
        break;
    }
    if (!is_fake)
      v_par_no_fake.push_back(p);
    else
      v_par_fake.push_back(p);
  }
  inside_polygons.push_back(to_polygon(v_par_fake));
  return v_par_no_fake;
}
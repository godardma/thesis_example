#include <codac>
#include <codac-capd.h>
#include "cgal_tools.h"

using namespace std;
using namespace codac2;

int main()
{
  // set_nb_threads(max_threads());
  
  ColorMap cmap = ColorMap::rainbow();
  VectorVar X(1);

  capd::IMap vectorField("var:x1,x2;fun:1-2*x1+1.5*sqr(x1)*x2,x1-1.5*sqr(x1)*x2;");

  double tf = 1.;
  double dt = tf;
  
  AnalyticFunction psi0 ({X},{0.05*X[0],0.05});
  OctaSym id ({1,2});
  OctaSym s ({-2,1});

  auto start_time = std::chrono::high_resolution_clock::now();

  auto peibos_output = PEIBOS(vectorField, tf, dt, psi0, {id,s,s*s,s.invert()}, 0.5, {0.95,0.05});

  auto m_v_par = reach_set(peibos_output);

  auto v_par = m_v_par[tf];

  Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_par);

  auto polygon_outer = to_codac(oriented_polygon(poly_with_holes.outer_boundary()));
  auto polygon_inner = to_codac(oriented_polygon(*poly_with_holes.holes_begin()));

  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start_time;
  printf("Computation time: %.4fs\n\n", elapsed.count());

  Figure2D output ("cora_inner",GraphicOutput::VIBES | GraphicOutput::IPE);
  output.set_axes(axis(0,{0.65,0.75}),axis(1,{0.485,0.585}));
  // output.set_axes(axis(0,{0.87,0.94}),axis(1,{0.74,0.81}));
  output.set_window_properties({100,100},{800,800});

  for (const auto& par: v_par)
    output.draw_parallelepiped(par, StyleProperties::boundary());

  output.draw_polygon(polygon_outer,StyleProperties(Color::blue(),"w:0.0002"));
  output.draw_polygon(polygon_inner,StyleProperties(Color::red(),"w:0.0002"));
}

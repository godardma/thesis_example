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

  capd::IMap vectorField("var:x1,x2;fun:1-2*x1+1.5*x1*x1*x2,x1-1.5*x1*x1*x2;");

  double tf = 1.0;
  double dt = 1.0;
  
  AnalyticFunction psi0 ({X},{0.125*X[0],0.125});
  OctaSym id ({1,2});
  OctaSym s ({-2,1});

  auto peibos_output = PEIBOS(vectorField, tf, dt, psi0, {id,s,s*s,s.invert()}, 0.25, {0.875,0.125}, true);

  auto m_v_par = reach_set(peibos_output);

  auto v_par = m_v_par[1.0];

  Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_par);

  auto polygon_outer = to_codac(oriented_polygon(poly_with_holes.outer_boundary()));
  auto polygon_inner = to_codac(oriented_polygon(*poly_with_holes.holes_begin()));

  Figure2D output ("cora_inner",GraphicOutput::VIBES | GraphicOutput::IPE);
  output.set_axes(axis(0,{0.6,0.8}),axis(1,{0.46,0.64}));
  output.set_window_properties({100,100},{800,800});

  // for (const auto& par: v_par)
  //   output.draw_parallelepiped(par, StyleProperties::boundary());

  output.draw_polygon(polygon_outer,Color::blue());
  output.draw_polygon(polygon_inner,Color::red());
}

#include <codac>
#include <codac-capd.h>

using namespace std;
using namespace codac2;

int main()
{
  set_nb_threads(max_threads());
  
  VectorVar X_3d(2);
  AnalyticFunction psi0_3d ({X_3d},{1,X_3d[0],X_3d[1]});

  capd::IMap vectorField_2d("var:x1,x2, x3; fun:cos(x3), sin(x3), 0;");

  OctaSym id_3d ({1,2,3});
  OctaSym s1 ({-2,1,3});
  OctaSym s2 ({3,2,-1});
  
  Figure3D figure3d ("Boat");
  figure3d.draw_axes(3.);

  Figure2D figure_3d_proj ("Boat projected", GraphicOutput::VIBES | GraphicOutput::IPE);
  figure_3d_proj.set_window_properties({25,25},{500,500});
  figure_3d_proj.set_axes({0,{-5,25}}, {1,{-20,20}});

  double tf = 20.;
  double dt = 5.;

  auto peibos_output = PEIBOS(vectorField_2d, tf, dt, psi0_3d, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, 0.125, true);  

  auto m_v_par = reach_set(peibos_output);

  for (auto& [t, v_par_2d] : m_v_par)
    for (const auto& par: v_par_2d)
      {
        figure3d.draw_parallelepiped(par, ColorMap::rainbow(0.5).color(t/tf));
        figure_3d_proj.draw_zonotope(par.proj({0,1}) , {Color::none(),ColorMap::rainbow(0.2).color(t/tf)});
      }
}

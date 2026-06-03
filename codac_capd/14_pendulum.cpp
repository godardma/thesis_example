#include <codac>
#include <codac-capd.h>

using namespace std;
using namespace codac2;

int main()
{
  set_nb_threads(max_threads());
  
  ColorMap cmap = ColorMap::rainbow();
  VectorVar X(1);

  capd::IMap vectorField("par:l,g;var:t,w;fun:w,-sin(t)*g/l - 0.5*w;");
 
  vectorField.setParameter("l",capd::interval(2.));
  vectorField.setParameter("g",capd::interval(10.));

  double tf = 10.0;
  double dt = 0.2;
  
  AnalyticFunction psi0 ({X},{0.1*sin(X[0]*PI/4),0.1*cos(X[0]*PI/4)});
  OctaSym id ({1,2});
  OctaSym s ({-2,1});

  auto peibos_output = PEIBOS(vectorField, tf, dt, psi0, {id,s,s*s,s.invert()}, 0.25, {-M_PI/4.,0.}, true);

  Figure2D output ("Pendulum",GraphicOutput::VIBES | GraphicOutput::IPE);
  output.set_axes(axis(0,{-1,0.8}),axis(1,{-1.2,1.8}));
  output.set_window_properties({800,100},{800,800});

  auto m_v_par = reach_set(peibos_output);

  for (auto& [t, v_par] : m_v_par)
    for (const auto& par: v_par)
      output.draw_parallelepiped(par, cmap.color(t/tf));

}
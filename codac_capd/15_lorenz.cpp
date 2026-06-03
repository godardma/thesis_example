#include <codac>
#include <codac-capd.h>

using namespace std;
using namespace codac2;

int main()
{
  set_nb_threads(max_threads());
  
  Figure3D figure ("Lorenz");
  figure.draw_axes(1.5);
  ColorMap cmap = ColorMap::rainbow(0.5);

  capd::IMap vectorField("par:sigma,rho,beta;var:x1,x2,x3;fun:10*(x2-x1),28*x1-x2-x1*x3,-2.6*x3+x1*x2;");
  vectorField.setParameter("sigma", 10.);
  vectorField.setParameter("rho", 28.);
  vectorField.setParameter("beta", 8/3);

  double tf = 0.1;
  double dt = 0.05;

  VectorVar X(2);
  AnalyticFunction psi0 ({X},{1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

  OctaSym id ({1,2,3});
  OctaSym s1 ({-2,1,3});
  OctaSym s2 ({3,2,-1});

  auto peibos_output = PEIBOS(vectorField, tf, dt, psi0, {id,s1,s1*s1,s1.invert(),s2,s2.invert()}, 0.1, true);

  auto m_v_par = reach_set(peibos_output);

  for (auto& [t, v_par] : m_v_par)
    for (const auto& par: v_par)
      figure.draw_parallelepiped(par, {cmap.color(t/tf),"t:"+to_string(t)});

}
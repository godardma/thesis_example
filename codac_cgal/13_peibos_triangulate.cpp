#include "triangulate_tools.h"

using namespace std;
using namespace codac2;

int main()
{
  set_nb_threads(max_threads());
  
  VectorVar y(2);
  AnalyticFunction f({y},{sqr(y[0])-sqr(y[1])+y[0],2*y[0]*y[1]+y[1]});

  VectorVar X(1);
  AnalyticFunction psi0 ({X},{sin(X[0]*PI/4.),cos(X[0]*PI/4.)});

  OctaSym id ({1,2});
  OctaSym s ({-2,1});

  auto v_par = PEIBOS(f, psi0, {id,s,s*s,s.invert()}, 0.125, true);

  IntervalVector y0({{-1.5,2.5},{-2,2}});

  triangulate_and_eliminate_fake(v_par, y0, "PEIBOS_triangulate");
}
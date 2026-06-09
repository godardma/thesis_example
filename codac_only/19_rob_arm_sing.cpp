// Author: Maël GODARD

#include<codac>

using namespace codac2;
using namespace std;

ColorMap rainbow_05(float alpha=1.)
{
  ColorMap cmap( Model::HSV );
  int i = 0;
  for(int h = 300 ; h > 0 ; h-=10)
  {
    cmap[i]=Color({(float)h,50.,100.,100.f*alpha},Model::HSV);
    i++;
  }
  return cmap;
}

int main()
{

  set_nb_threads(max_threads());

  VectorVar y(3);
  double l1=0.45;
  double l2=0.3;
  double l3=0.3;
  AnalyticFunction f({y},{cos(y[0])*(l1+l2*cos(y[1])+l3*cos(y[1]+y[2])),sin(y[0])*(l1+l2*cos(y[1])+l3*cos(y[1]+y[2])),l2*sin(y[1])+l3*sin(y[1]+y[2])});

  
  VectorVar X(2);
  AnalyticFunction psi0 ({X},{PI/2.,X[0]*PI/2.,X[1]*PI/2.});

  OctaSym id_3d ({1,2,3});
  OctaSym s1 ({-2,1,3});
  OctaSym s2 ({3,2,-1});

  double epsilon = 0.05;

  vector<OctaSym> Sigma({id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()});

  Figure3D output ("rob_arm_sing");

  for (int i = 0; i < Sigma.size(); i++)
  {
    auto v_par_1 = PEIBOS(f, psi0, {Sigma[i]}, epsilon);
    for (auto p : v_par_1)
      output.draw_parallelepiped(p,rainbow_05(0.5).color((double)(i)/(double)(Sigma.size()-1)));
  }

  IntervalVector X0 ({{-1,1},{-1,1}});
  OctaSym symmetry ({1,2,3});

  AnalyticFunction psi ({X},{X[0]*PI/2.,X[1]*PI/2.,0.});

  auto v_par_2 = PEIBOS(f, psi, {symmetry}, epsilon, {0.0,0.0,0.0});
  for (auto p : v_par_2)
    output.draw_parallelepiped(p,Color::blue(0.7));

  output.draw_axes(0.5,{0.2,-1.2,-0.8});
}

  

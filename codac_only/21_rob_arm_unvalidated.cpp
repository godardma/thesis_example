// Author: Maël GODARD

#include <codac>
#include <codac-unsupported.h>

#include <chrono>

using namespace std;
using namespace codac2;

bool verify_criteria (Parallelepiped& p)
{
  CtcPolar c;
  IntervalVector bbox_polar = p.box();
  Interval r(1.05,+oo), theta(0,2.0*M_PI);
  c.contract(bbox_polar[1],bbox_polar[2], r, theta);

  IntervalVector bbox_box = p.box();
  Interval obstacle (1.055,+oo);
  CtcWrapper ctc(obstacle);
  ctc.contract(bbox_box[0]);

  return bbox_polar.is_empty() && bbox_box.is_empty();
}

template <typename T>
void PEIBOS_recursive (AnalyticFunction<T>& f, AnalyticFunction<T>& psi_0, OctaSym symmetry , IntervalVector X, double epsilon, Figure3D& figure_3d)
{

  VectorVar x(2);
  AnalyticFunction g ({x},f(symmetry(psi_0(x))));
  AnalyticFunction g_id ({x},symmetry(psi_0(x)));
  if (epsilon < 0.01)
  {
    auto parallel = g.parallelepiped_eval(X);
    figure_3d.draw_parallelepiped(parallel, StyleProperties(Color::red(0.8),"workspace"));

    AnalyticFunction g_id ({x},symmetry(psi_0(x)));
    auto parallel_id = g_id.parallelepiped_eval(X);
    figure_3d.draw_parallelepiped(parallel_id, StyleProperties(Color::red(),"atlas"));
  }

  else
  {
    std::vector<IntervalVector> boxes;
    double true_eps = split(X, epsilon, boxes);
    for (auto& X0 : boxes)
    {
        auto parallel = g.parallelepiped_eval(X0);
        if (verify_criteria(parallel))
        {
          figure_3d.draw_parallelepiped(parallel, StyleProperties(Color::green(0.5),"workspace"));
          auto parallel_id = g_id.parallelepiped_eval(X);
          figure_3d.draw_parallelepiped(parallel_id, StyleProperties(Color::green(0.2),"atlas"));
        }
        else
          PEIBOS_recursive(f, psi_0, symmetry, X0, epsilon/2., figure_3d);
    }
  }
}

template <typename T>
void PEIBOS_adaptative(AnalyticFunction<T>& f, AnalyticFunction<T>& psi_0,const vector<OctaSym>& symmetries , double epsilon, Figure3D& figure_3d)
{
  double epsilon_init = epsilon;

  for (int i = 0; i < symmetries.size(); i++)
    {
      IntervalVector X = IntervalVector({{-1,1},{-1,1}});
      PEIBOS_recursive(f, psi_0, symmetries[i], X, epsilon, figure_3d);
    }
}

int main()
{  
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

  double epsilon = 1.;

  Figure3D output ("rob_arm_nonvalidated");
  output.draw_axes(0.5,{0,-1.2,-0.8});

  vector<Parallelepiped> v_par ;

  IntervalVector X0 ({{-1,1},{-1,1}});
  OctaSym symmetry ({1,2,3});

  AnalyticFunction psi ({X},{X[0]*PI/2.,X[1]*PI/2.,0.});

  // start the timer
  auto start = chrono::high_resolution_clock::now();

  PEIBOS_adaptative(f, psi0, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, epsilon, output);
  PEIBOS_recursive(f, psi, {symmetry}, X0, epsilon, output);

    // stop the timer
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
  cout << "Execution time: " << duration.count()/1000. << " seconds" << endl;

}

  

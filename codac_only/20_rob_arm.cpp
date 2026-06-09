// Author: Maël GODARD

#include <codac>
#include <codac-unsupported.h>

#include <chrono>

using namespace std;
using namespace codac2;

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

bool verify_criteria (Parallelepiped& p)
{
  CtcPolar c;
  IntervalVector bbox_polar = p.box();
  Interval r(1.2,+oo), theta(0,2.0*M_PI);
  c.contract(bbox_polar[1],bbox_polar[2], r, theta);

  IntervalVector bbox_box = p.box();
  Interval obstacle (1.055,+oo);
  CtcWrapper ctc(obstacle);
  ctc.contract(bbox_box[0]);

  return bbox_polar.is_empty() && bbox_box.is_empty();
}

template <typename T>
void PEIBOS_recursive (AnalyticFunction<T>& g , IntervalVector X, double epsilon, Figure3D& figure_3d, vector<Parallelepiped>& output)
{
  VectorVar x(3);
  AnalyticFunction f_id({x},x);
  if (epsilon < 0.01)
    cout << "epsilon limit reached" << endl;
  else
  {
    std::vector<IntervalVector> boxes;
    double true_eps = split(X, epsilon, boxes);
    for (auto& X0 : boxes)
    {
        auto parallel = g.parallelepiped_eval(X0);
        if (verify_criteria(parallel))
        {
          output.push_back(Parallelepiped(parallel.c,parallel.A));
          figure_3d.draw_parallelepiped(parallel, rainbow_05(0.5).color(epsilon/0.5));
        }
        else
          PEIBOS_recursive(g, X0, epsilon/2., figure_3d, output);
      }
  }
}

template <typename T>
void PEIBOS_adaptative(AnalyticFunction<T>& f, AnalyticFunction<T>& psi_0,const vector<OctaSym>& symmetries , double epsilon, Figure3D& figure_3d, vector<Parallelepiped>& output)
{
  double epsilon_init = epsilon;
  for (int i = 0; i < symmetries.size(); i++)
    {
      VectorVar x(2);
      AnalyticFunction g ({x},f(symmetries[i](psi_0(x))));
      IntervalVector X = IntervalVector({{-1,1},{-1,1}});
      PEIBOS_recursive(g, X, epsilon, figure_3d,output);
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

  Figure3D output ("rob_arm");

  output.draw_axes(0.5,{0,-1.2,-0.8});

  vector<Parallelepiped> v_par ;

  IntervalVector X0 ({{-1,1},{-1,1}});
  OctaSym symmetry ({1,2,3});

  AnalyticFunction psi ({X},{X[0]*PI/2.,X[1]*PI/2.,0.});

  // Adaptive /////////////////////////////

  // start the timer
  auto start = chrono::high_resolution_clock::now();
  PEIBOS_adaptative(f, psi0, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, epsilon, output, v_par);

  VectorVar x(2);
  AnalyticFunction g_sing ({x},f(symmetry(psi(x))));
  PEIBOS_recursive(g_sing, X0, epsilon, output, v_par);

    // stop the timer
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
  cout << "Execution time for adaptive PEIBOS: " << duration.count()/1000. << " seconds" << endl;

  ////////////////////////////////////////

  // Classical ///////////////////////////

  // start the timer
  start = chrono::high_resolution_clock::now();

  auto v_par_1 = PEIBOS(f, psi0, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, 0.03125, {0.0,0.0,0.0});
  auto v_par_2 = PEIBOS(f, psi, {symmetry}, 0.03125, {0.0,0.0,0.0});


  // stop the timer
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
  cout << "Execution time: for classical PEIBOS " << duration.count()/1000. << " seconds" << endl;

  v_par_1.insert(v_par_1.end(), v_par_2.begin(), v_par_2.end());

  ////////////////////////////////////////

  IntervalVector X_limit ({{1.055,1.5},{-2,2}});


  Figure2D output_xy ("rob_arm_proj_xy_adapt", GraphicOutput::VIBES | GraphicOutput::IPE);
  output_xy.set_window_properties({50,50},{500,500});
  output_xy.set_axes(axis(0,{-0.5,1.5}),axis(1,{-2,2}));
  output_xy.draw_box(X_limit,{Color({255,102,0,255}),Color({255,102,0,130})});


  Figure2D output_yz ("rob_arm_proj_yz_adapt", GraphicOutput::VIBES | GraphicOutput::IPE);
  output_yz.set_window_properties({650,50},{500,500});
  output_yz.set_axes(axis(0,{-2,2}),axis(1,{-2,2}));
  output_yz.draw_ring({0.,0.},{1.2,2},{Color({255,102,0,255}),Color({255,102,0,130})});

  Figure2D output_xy_fine ("rob_arm_proj_xy_fine", GraphicOutput::VIBES | GraphicOutput::IPE);
  output_xy_fine.set_window_properties({50,650},{500,500});
  output_xy_fine.set_axes(axis(0,{-0.5,1.5}),axis(1,{-2,2}));
  output_xy_fine.draw_box(X_limit,{Color({255,102,0,255}),Color({255,102,0,130})});


  Figure2D output_yz_fine ("rob_arm_proj_yz_fine", GraphicOutput::VIBES | GraphicOutput::IPE);
  output_yz_fine.set_window_properties({650,650},{500,500});
  output_yz_fine.set_axes(axis(0,{-2,2}),axis(1,{-2,2}));
  output_yz_fine.draw_ring({0.,0.},{1.2,2},{Color({255,102,0,255}),Color({255,102,0,130})});


  for (auto p : v_par)
  {
    output_xy.draw_zonotope(p.proj({0,1}), {Color::black(), Color::green(0.5)});
    output_yz.draw_zonotope(p.proj({1,2}), {Color::black(), Color::green(0.5)});
  }

  for (auto p : v_par_1)
  {
    output_xy_fine.draw_zonotope(p.proj({0,1}), {Color::black(), Color::green(0.5)});
    output_yz_fine.draw_zonotope(p.proj({1,2}), {Color::black(), Color::green(0.5)});
  }
}

  

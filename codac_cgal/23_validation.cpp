#include "triangulate_tools.h"
#include <iostream>
#include <chrono>

using namespace std;
using namespace codac2;

vector<Parallelepiped> v_par_ad_all;
vector<IntervalVector> v_boxes_in;
double r;

bool verify_criteria (const Parallelepiped& p)
{
  VectorVar X(2);
  auto decomp = inverse_enclosure(p.A)*(X-p.c);

  AnalyticFunction f_polar({X},{decomp[0],decomp[1], sqrt(sqr(X[0]-0.5)+sqr(X[1]))});
  IntervalVector value_polar ({{-1,1},{-1,1},{1.75,oo}});
  CtcInverse ctc_inverse_polar(f_polar,value_polar);
  
  IntervalVector X1 = p.box();
  ctc_inverse_polar.contract(X1);
  return X1.is_empty();
}

template <typename T>
void PEIBOS_recursive (AnalyticFunction<T>& g , IntervalVector X, double epsilon, vector<Parallelepiped>& output,vector<Parallelepiped>& output_invalid)
{
  if (epsilon < 0.0625)
  {
    auto parallel = g.parallelepiped_eval(X);
    output_invalid.push_back(parallel);
    v_par_ad_all.push_back(parallel);
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
          v_par_ad_all.push_back(parallel);
          output.push_back(parallel);
        }
        else
          PEIBOS_recursive(g, X0, epsilon/2., output,output_invalid);
      }
  }
}

template <typename T>
void PEIBOS_adaptative(AnalyticFunction<T>& f, AnalyticFunction<T>& psi_0,const vector<OctaSym>& symmetries , double epsilon, vector<Parallelepiped>& output,vector<Parallelepiped>& output_invalid)
{
  double epsilon_init = epsilon;
  for (int i = 0; i < symmetries.size(); i++)
    {
      VectorVar x(1);
      AnalyticFunction g ({x},f(symmetries[i](psi_0(x))));
      IntervalVector X = IntervalVector({{-1,1}});
      PEIBOS_recursive(g, X, epsilon,output,output_invalid);
    }
}

BoolInterval test_inside(const Vector& v)
{    
  return BoolInterval::FALSE;
}

int main()
{
  r = 1.75;

  VectorVar x (2);
  AnalyticFunction c ({x},sqrt(sqr(x[0])+sqr(x[1]))-1);

  AnalyticFunction f ({x},{sqr(x[0])-sqr(x[1])+x[0],2*x[0]*x[1]+x[1]});
  AnalyticFunction f_circ ({x},{r*x[0]+0.5,r*x[1]});

  AnalyticFunction g ({x},sqrt(sqr(x[0]-0.5)+sqr(x[1]))-r);

  AnalyticFunction h ({x},g(f(x)));

  AnalyticFunction ch ({x},{c(x),h(x)});

  VectorVar X(1);
  AnalyticFunction psi0 ({X},{sin(X[0]*PI/4.),cos(X[0]*PI/4.)});

  OctaSym id ({1,2});
  OctaSym s ({-2,1});
  vector<OctaSym> symmetries({id,s,s*s,s.invert()});

  // PEIBOS
  auto v_par = PEIBOS(f, psi0, symmetries, 0.0625, true);

  // Inversion
  IntervalVector X0 = IntervalVector::constant(2,{-1.5,1.5});
  IntervalVector Y0 ({{-1.5,2.5},{-2.,2.}});

  Figure2D fig_inversion_in_out("inversion_in_out", GraphicOutput::VIBES | GraphicOutput::IPE);

  fig_inversion_in_out.set_axes(Y0);
  fig_inversion_in_out.set_window_properties({1150, 50}, {500, 500});
  fig_inversion_in_out.draw_circle({0.5,0},1.75,StyleProperties(Color::blue(),"z:-5"));

  auto start_time = std::chrono::high_resolution_clock::now();

  SepInverse sep_inv_out (ch,IntervalVector({{-oo,0},{0,oo}}));

  auto p_out = pave (X0,sep_inv_out,0.01);
  auto v_cs = p_out.connected_subsets(PavingInOut::bound);
  CtcUnion ctc_union(2);
  vector<vector<IntervalVector>> y_out;
  int i = 0;
  for (const auto & cs:v_cs)
  {
    auto x_out = cs.boxes();
    vector<IntervalVector> y_out_i;
    for (const auto& x : x_out)
    {
      auto y = f.eval(x);
      y_out_i.push_back(y);
      ctc_union|=CtcWrapper<IntervalVector>(f.eval(x));
    }
    y_out.push_back(y_out_i);
  }
  
  for (const auto& y_out_i : y_out)
    triangulate(y_out_i, fig_inversion_in_out);

  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start_time;
  printf("Computation time for inversion in and out: %.4fs\n\n", elapsed.count());

  // Ad-PEIBOS
  vector<Parallelepiped> v_par_ad;
  vector<Parallelepiped> v_par_ad_invalid;

  start_time = std::chrono::high_resolution_clock::now();

  PEIBOS_adaptative(f,psi0,symmetries,2.,v_par_ad, v_par_ad_invalid);
  auto v_par_circle = PEIBOS(f_circ, psi0, symmetries, 0.2);

  triangulate_with_diff(v_par_ad_all,v_par_circle,Y0,"validation");

  elapsed = std::chrono::high_resolution_clock::now() - start_time;
  printf("Computation time for Ad-PEIBOS: %.4fs\n\n", elapsed.count());

  Figure2D fig_inversion_x("inversion_x", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_inversion_x.set_window_properties({50,50},{500,500});
  fig_inversion_x.set_axes(X0);

  Figure2D fig_reference ("reference", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_reference.set_window_properties({50,50},{500,500});
  fig_reference.set_axes(Y0);

  Figure2D fig_inversion ("inversion", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_inversion.set_window_properties({600,50},{500,500});
  fig_inversion.set_axes(Y0);

  Figure2D fig_PEIBOS ("AD-PEIBOS", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_PEIBOS.set_window_properties({50,600},{500,500});
  fig_PEIBOS.set_axes(Y0);

  for (const auto& y_out_i : y_out)
    for (const auto& y : y_out_i)
      fig_inversion.draw_box(y,StyleProperties::boundary());

  for (const auto& p : v_par_ad)
    fig_PEIBOS.draw_parallelepiped(p,Color::green());

  for (const auto& p : v_par_ad_invalid)
    fig_PEIBOS.draw_parallelepiped(p,Color::red());

  for (const auto& p : v_par)
    fig_reference.draw_parallelepiped(p,Color::black());
  fig_inversion_x.draw_paving(p_out);

  fig_reference.draw_circle({0.5,0},r,Color::blue());
  fig_inversion_x.draw_circle({0,0},1.0,Color::black());
  fig_inversion.draw_circle({0.5,0},r,Color::blue());
  fig_PEIBOS.draw_circle({0.5,0},r,Color::blue());


}
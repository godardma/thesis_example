
#include <codac>
#include <iostream>
#include <chrono>

using namespace std;
using namespace codac2;

vector<Parallelepiped> v_par_ad_all;
vector<IntervalVector> v_boxes_in;

bool verify_criteria (const Parallelepiped& p)
{
  // CtcDist c;
  // IntervalVector bbox_dist = p.box();
  // Interval r(1.8,+oo);
  // IntervalVector x_dist ({{0.5,0.5}, {0,0}, bbox_dist[0], bbox_dist[1] , r});
  // c.contract(x_dist);

  // return x_dist.is_empty();

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
    // cout << "epsilon limit reached" << endl;
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

// bool verify_criteria (const IntervalVector& bbox_dist)
// {
//   CtcDist c;
//   Interval r(1.8,+oo);
//   IntervalVector x_dist ({{0.5,0.5}, {0,0}, bbox_dist[0], bbox_dist[1] , r});
//   c.contract(x_dist);

//   return x_dist.is_empty();
// }


BoolInterval test_inside(const Vector& v)
{    

  // for (const auto& p_in : v_points_in)
  // {
  //   if (sqr(v[0]-p_in[0])+sqr(v[1]-p_in[1]) < 0.01)
  //     return BoolInterval::TRUE;
  // }

  for (const auto& p_in : v_boxes_in)
  {
    if (p_in.contains(v))
      return BoolInterval::TRUE;
  }

  return BoolInterval::FALSE;
}

int main()
{
  double r = 1.75;

  VectorVar x (2);
  AnalyticFunction c ({x},sqrt(sqr(x[0])+sqr(x[1]))-1);

  AnalyticFunction f ({x},{sqr(x[0])-sqr(x[1])+x[0],2*x[0]*x[1]+x[1]});
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

  // Ad-PEIBOS
  vector<Parallelepiped> v_par_ad;
  vector<Parallelepiped> v_par_ad_invalid;
  auto start_time = std::chrono::high_resolution_clock::now();
  PEIBOS_adaptative(f,psi0,symmetries,2.,v_par_ad, v_par_ad_invalid);
  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start_time;
  printf("Computation time for Ad-PEIBOS: %.4fs\n\n", elapsed.count());

  // Inversion
  IntervalVector X0 = IntervalVector::constant(2,{-1.5,1.5});
  IntervalVector Y0 ({{-2,3},{-2.5,2.5}});

  start_time = std::chrono::high_resolution_clock::now();

  SepInverse sep_inv_out (ch,IntervalVector({{-oo,0},{0,oo}}));

  auto p_out = pave (X0,sep_inv_out,0.01);
  auto x_out = p_out.boxes(PavingInOut::inner);
  vector<IntervalVector> y_out;
  for (const auto& x : x_out)
    y_out.push_back(f.eval(x));

  elapsed = std::chrono::high_resolution_clock::now() - start_time;
  printf("Computation time for inversion: %.4fs\n\n", elapsed.count());

  // Paving circe
  IntervalVector polar_constraint({{0,r},{0,2*PI}});
  SepWrapper<IntervalVector> sep_wrapper (polar_constraint);
  SepPolarCart sep_cart_polar (sep_wrapper);
  SepTransform sep_transform(sep_cart_polar,
      AnalyticFunction({x}, x+Vector({0.5,0})),
      AnalyticFunction({x}, x-Vector({0.5,0})));
  auto p_circ = pave(Y0,sep_transform,0.01);

  // Paving AD-PEIBOS output
  CtcUnion ctc_union(2);
  for (const auto& p : v_par_ad_all)
    ctc_union|=CtcWrapper<Parallelepiped>(p);

  SepInverse sep_unit_circle (c,Interval({-oo,0}));
  auto p_unit_circle = pave(X0,sep_unit_circle,0.05);
  auto v_b_inner = p_unit_circle.boxes(PavingInOut::inner);
  for (const auto& b : v_b_inner)
    v_boxes_in.push_back(f.eval(b));

  SepCtcBoundary sep_ctc_boundary(ctc_union,test_inside);

  Figure2D fig_reference ("reference", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_reference.set_window_properties({50,50},{500,500});
  fig_reference.set_axes(Y0);

  Figure2D fig_inversion ("inversion", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_inversion.set_window_properties({600,50},{500,500});
  fig_inversion.set_axes(Y0);

  Figure2D fig_PEIBOS ("AD-PEIBOS", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_PEIBOS.set_window_properties({1150,50},{500,500});
  fig_PEIBOS.set_axes(Y0);

  Figure2D fig_paving_circle ("cirlce_paving", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_paving_circle.set_window_properties({50,600},{500,500});
  fig_paving_circle.set_axes(Y0);

  Figure2D fig_paving_conform ("conform_paving", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_paving_conform.set_window_properties({600,600},{500,500});
  fig_paving_conform.set_axes(Y0);

  Figure2D fig_inversion_x("inversion_x", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_inversion_x.set_window_properties({1150,600},{500,500});
  fig_inversion_x.set_axes(X0);



  // fig_paving_conform.pave(Y0,sep_ctc_boundary,0.1);
  auto p_boundary = pave(Y0,sep_ctc_boundary,0.4);
  fig_paving_conform.draw_paving(p_boundary);

  for (const auto& p_in : v_boxes_in)
    fig_paving_conform.draw_box(p_in, Color::red());

  for (const auto& y : y_out)
    fig_inversion.draw_box(y,StyleProperties::boundary());

  for (const auto& p : v_par_ad)
    fig_PEIBOS.draw_parallelepiped(p,Color::green());

  for (const auto& p : v_par_ad_invalid)
    fig_PEIBOS.draw_parallelepiped(p,Color::red());

  for (const auto& p : v_par)
    fig_reference.draw_parallelepiped(p,Color::black());

  fig_paving_circle.draw_paving(p_circ);
  fig_reference.draw_circle({0.5,0},r,Color::blue());
  fig_inversion.draw_circle({0.5,0},r,Color::blue());
  fig_PEIBOS.draw_circle({0.5,0},r,Color::blue());
  // fig_IEIBOS.draw_circle({0.5,0},r,Color::blue());

  fig_inversion_x.draw_paving(p_out);

}
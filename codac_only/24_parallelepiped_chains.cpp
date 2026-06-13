#include <codac>
#include <vector>

using namespace codac2;
using namespace std;

class ParallelepipedLink
{
  public:
    AnalyticFunction<VectorType> psi;
    IntervalVector x;
    bool is_intersected = false;
    ParallelepipedLink(const AnalyticFunction<VectorType>& psi, const IntervalVector& x) : psi(psi), x(x) {}

    bool is_neighbor(const ParallelepipedLink& other) const
    {
      return x.intersects(other.x);
    }
};

class ParallelepipedChain : public vector<ParallelepipedLink>
{
  public:
    bool is_neighbor(const ParallelepipedLink& link)
    {
      for (const auto& l : *this)
        if (l.is_neighbor(link))
          return true;
      return false;
    }

    IntervalVector bounding_box() const
    {
      IntervalVector box = (*this)[0].x;
      for (const auto& link : *this)
        box |= link.x;
      return box;
    }
};

bool injectivity_criterion(const AnalyticFunction<VectorType>& f, const IntervalVector& x)
{
  IntervalMatrix Jf = f.diff(x);
  IntvFullPivLU lu(Jf);
  return !lu.determinant().contains(0);
}

bool injectivity_criterion (const AnalyticFunction<VectorType>& f, const ParallelepipedChain& chain, const ParallelepipedLink& link)
{
  return injectivity_criterion(f,link.psi.eval(chain.bounding_box() | link.x));
}

vector<ParallelepipedChain> BC_decomposition(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi, const vector<IntervalVector>& cover)
{
  auto L_W = cover;
  vector<ParallelepipedChain> L_BC;
  IntervalVector x = L_W.back();
  L_W.pop_back();
  while (!L_W.empty())
  {
    ParallelepipedChain chain;
    chain.push_back(ParallelepipedLink(psi, x));
    while (!L_W.empty())
    {
      x = L_W.back();
      L_W.pop_back();
      ParallelepipedLink link (psi,x);
      if (injectivity_criterion(f,chain,link) && chain.is_neighbor(link))
        chain.push_back(link);
      else
        break;
    }
    L_BC.push_back(chain);
  }
  return L_BC;
}

bool find_intersection(const AnalyticFunction<VectorType>& f, vector<ParallelepipedChain>& L_BC)
{
  bool found_intersection = false;
  // for (size_t i = 0; i < L_BC.size(); ++i)
  // {
  //   auto& chain1 = L_BC[i];
  //   VectorVar X(chain1[0].x.size());
  //   AnalyticFunction g1 ({X},f(chain1[0].psi(X)));
  //   for (size_t j = i+1; j < L_BC.size(); ++j)
  //   {
  //     auto& chain2 = L_BC[j];
  //     AnalyticFunction g2 ({X},f(chain2[0].psi(X)));
  //     for (auto& link1 : chain1)
  //     {
  //       for (auto& link2 : chain2)
  //       {
  //         auto z1 = g1.eval(link1.x);
  //         auto z2 = g2.eval(link2.x);
  //         if (z1.intersects(z2))
  //         {
  //           auto y1 = link1.psi.eval(link1.x);
  //           auto y2 = link2.psi.eval(link2.x);
  //           if (!injectivity_criterion(f, y1|y2))
  //           {
  //             found_intersection = true;
  //             link1.is_intersected = true;
  //             link2.is_intersected = true;
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  return found_intersection;
}

int main()
{
  // Graphical setup

  ColorMap cmap = ColorMap::rainbow(0.5);


  Figure3D fig_torus("torus_atlas");
  fig_torus.draw_axes(0.7);
  // Constructing covers

  IntervalVector V1 ({{-1,1},{-1,1}});
  IntervalVector V2 ({{-1,1},{-1,1}});
  IntervalVector V3 ({{-1,1},{-1,1}});
  IntervalVector V4 ({{-1,1},{-1,1}});

  vector<IntervalVector> C1;
  vector<IntervalVector> C2;
  vector<IntervalVector> C3;
  vector<IntervalVector> C4;

  double epsilon = 0.1;

  split(V1, epsilon, C1);
  split(V2, epsilon, C2);
  split(V3, epsilon, C3);
  split(V4, epsilon, C4);

  cout << "Number of boxes in cover: " << C1.size() << endl;

  // Box Chain decomposition

  double r = 1.;
  double R = 1.8;

  VectorVar X(2);
  // gives a quarter of the torus
  AnalyticFunction psi_1 ({X},{(R + r*sin(X[0]*PI/2))*cos(X[1]*PI/2),
                          (R + r*sin(X[0]*PI/2))*sin(X[1]*PI/2),
                          r*cos(X[0]*PI/2)});
  AnalyticFunction psi_2 ({X},{(R + r*sin(PI+X[0]*PI/2))*cos(X[1]*PI/2),
                          (R + r*sin(PI+X[0]*PI/2))*sin(X[1]*PI/2),
                          r*cos(PI+X[0]*PI/2)});

  AnalyticFunction psi_3 ({X},{(R + r*sin(X[0]*PI/2))*cos(PI+X[1]*PI/2),
                          (R + r*sin(X[0]*PI/2))*sin(PI+X[1]*PI/2),
                          r*cos(X[0]*PI/2)});

  AnalyticFunction psi_4 ({X},{(R + r*sin(PI+X[0]*PI/2))*cos(PI+X[1]*PI/2),
                          (R + r*sin(PI+X[0]*PI/2))*sin(PI+X[1]*PI/2),
                          r*cos(PI+X[0]*PI/2)});

  ScalarVar t;
  VectorVar x(2);

  AnalyticFunction traj({t},vec(-9.6* sqr(t)  + 11.4*(0.5* pow(t,3) - t) +30,-9.6*(0.5* pow(t,3) - t) - 11.4* sqr(t) + 30));
  AnalyticFunction d_traj({t},vec(-19.2* t + 17.1* sqr(t) - 11.4,-14.4* sqr(t) +9.6 -22.8* t));

  AnalyticFunction f_loop({x},vec(-9.6* sqr(x[1])  + 11.4*(0.5* pow(x[1],3) - x[1]) +30 + x[0]*sin(atan2(-14.4* sqr(x[1]) +9.6 -22.8* x[1],-19.2* x[1] + 17.1* sqr(x[1]) - 11.4)),-9.6*(0.5* pow(x[1],3) - x[1]) - 11.4* sqr(x[1]) + 30-x[0]*cos(atan2(-14.4* sqr(x[1]) +9.6 -22.8* x[1],-19.2* x[1] + 17.1* sqr(x[1]) - 11.4))));

  // for (const auto&box : C1)
  //   fig_torus.draw_parallelepiped(psi_1.parallelepiped_eval(box), Color::red(0.5));

  // for (const auto&box : C2)
  //   fig_torus.draw_parallelepiped(psi_2.parallelepiped_eval(box), Color::green(0.5));

  // for (const auto&box : C3)
  //   fig_torus.draw_parallelepiped(psi_3.parallelepiped_eval(box), Color::blue(0.5));

  // for (const auto&box : C4)
  //   fig_torus.draw_parallelepiped(psi_4.parallelepiped_eval(box), Color::orange(0.5));

  // VectorVar Y(2);
  // AnalyticFunction f ({Y},{sqr(Y[0])-sqr(Y[1])+Y[0],2*Y[0]*Y[1]+Y[1]});

  VectorVar Y (3);
  AnalyticFunction f ({Y},{Y[0],f_loop(Y.subvector(1,2))[0],f_loop(Y.subvector(1,2))[1]});

  AnalyticFunction g1 ({X},f(psi_1(X)));
  AnalyticFunction g2 ({X},f(psi_2(X)));
  AnalyticFunction g3 ({X},f(psi_3(X)));
  AnalyticFunction g4 ({X},f(psi_4(X)));
  
    for (const auto&box : C1)
    fig_torus.draw_parallelepiped(g1.parallelepiped_eval(box), Color::red(0.5));

  for (const auto&box : C2)
    fig_torus.draw_parallelepiped(g2.parallelepiped_eval(box), Color::green(0.5));

  for (const auto&box : C3)
    fig_torus.draw_parallelepiped(g3.parallelepiped_eval(box), Color::blue(0.5));

  for (const auto&box : C4)
    fig_torus.draw_parallelepiped(g4.parallelepiped_eval(box), Color::orange(0.5));

  // auto L_BC1 = BC_decomposition(f, psi_1, C1);
  // auto L_BC2 = BC_decomposition(f, psi_2, C2);

  // auto L_BC = L_BC1;
  // for (const auto& chain : L_BC2)
  //   L_BC.push_back(chain);

  // cout << "Number of box chains: " << L_BC.size() << endl;

  // for (size_t i = 0; i < L_BC.size(); ++i)
  // {
  //   const auto& chain = L_BC[i];
  //   for (const auto& link : chain)
  //   {
  //     AnalyticFunction gi ({X},f(link.psi(X)));
  //     double color_value = (double)i/(double)(L_BC.size()-1);
  //     fig_bc_init.draw_box(link.psi.eval(link.x), {cmap_edge.color(color_value), cmap_fill.color(color_value)});
  //     fig_bc_image.draw_box(gi.eval(link.x), {cmap_edge.color(color_value), cmap_fill.color(color_value)});
  //   }
  // }

  // // Detecting self-intersections

  // bool intersection_found = find_intersection(f, L_BC);
  // cout << "intersection found ? "<< (intersection_found?"true":"false") << endl;

  // vector<ParallelepipedLink> unintersected_links;
  // vector<ParallelepipedLink> intersected_links;
  // vector<ParallelepipedLink> L_BL;

  // for (const auto& chain : L_BC)
  // {
  //   for (const auto& link : chain)
  //   {
  //     L_BL.push_back(link);
  //     AnalyticFunction gi ({X},f(link.psi(X)));
  //     if (link.is_intersected)
  //     {
  //       intersected_links.push_back(link);
  //       fig_box_chains_new.draw_box(gi.eval(link.x), {Color::red(), Color::red(0.5)});
  //       fig_intersect.draw_box(gi.eval(link.x), {Color::red(), Color::red(0.5)});
  //       fig_fake_boundary_removed.draw_box(gi.eval(link.x), {Color::red(), Color::red(0.5)});

  //     }
  //     else
  //     {
  //       unintersected_links.push_back(link);
  //       fig_intersect.draw_box(gi.eval(link.x), {Color::black(), Color::black(0.5)});
  //     }
  //   }
  // }
}
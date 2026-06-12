#include <codac>
#include <vector>

using namespace codac2;
using namespace std;

class BoxLink
{
  public:
    AnalyticFunction<VectorType> psi;
    IntervalVector x;
    BoxLink(const AnalyticFunction<VectorType>& psi, const IntervalVector& x) : psi(psi), x(x) {}

    bool is_neighbor(const BoxLink& other) const
    {
      return x.intersects(other.x);
    }
};

class BoxChain : public vector<BoxLink>
{
  public:
    bool is_neighbor(const BoxLink& link)
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

bool injectivity_criterion (const AnalyticFunction<VectorType>& f, const BoxChain& chain, const BoxLink& link)
{
  IntervalVector new_bounding_box = chain.bounding_box() | link.x;
  IntervalMatrix Jf = f.diff(link.psi.eval(new_bounding_box));
  IntvFullPivLU lu(Jf);
  return !lu.determinant().contains(0);
}

vector<BoxChain> BC_decomposition(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi, const vector<IntervalVector>& cover)
{
  auto L_W = cover;
  vector<BoxChain> L_BC;
  IntervalVector x = L_W.back();
  L_W.pop_back();
  while (!L_W.empty())
  {
    BoxChain chain;
    chain.push_back(BoxLink(psi, x));
    while (!L_W.empty())
    {
      x = L_W.back();
      L_W.pop_back();
      BoxLink link (psi,x);
      if (injectivity_criterion(f,chain,link) && chain.is_neighbor(link))
        chain.push_back(link);
      else
        break;
    }
    L_BC.push_back(chain);
  }
  return L_BC;
}

bool find_intersection (const vector<BoxChain>& L_BC, const AnalyticFunction<VectorType>& f)
{
  return true;
}

int main()
{
  IntervalVector V1 ({{-1,1}});
  IntervalVector V2 ({{-1,1}});

  vector<IntervalVector> C1;
  vector<IntervalVector> C2;

  double epsilon = 0.1;

  split(V1, epsilon, C1);
  split(V2, epsilon, C2);

  cout << "Number of boxes in cover: " << C1.size() << endl;

  VectorVar X(1);
  AnalyticFunction psi_1 ({X},{cos(X[0]*PI/2),sin(X[0]*PI/2)});
  AnalyticFunction psi_2 ({X},{cos(PI+X[0]*PI/2),sin(PI+X[0]*PI/2)});

  VectorVar Y(2);
  AnalyticFunction f ({Y},{sqr(Y[0])-sqr(Y[1])+Y[0],2*Y[0]*Y[1]+Y[1]});

  AnalyticFunction g1 ({X},f(psi_1(X)));
  AnalyticFunction g2 ({X},f(psi_2(X)));

  auto L_BC1 = BC_decomposition(f, psi_1, C1);
  auto L_BC2 = BC_decomposition(f, psi_2, C2);

  auto L_BC = L_BC1;
  for (const auto& chain : L_BC2)
    L_BC.push_back(chain);

  cout << "Number of box chains: " << L_BC.size() << endl;

  ColorMap cmap_edge = ColorMap::rainbow();
  ColorMap cmap_fill = ColorMap::rainbow(0.5);

  IntervalVector axis_limit_initial ({{-1.5,1.5},{-1.5,1.5}});
  IntervalVector axis_limit_image ({{-2.25,2.75},{-2.5,2.5}});

  Figure2D fig_bc_init("box_chains_init", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_bc_init.set_window_properties({50,50},{500,500});
  fig_bc_init.set_axes(axis_limit_initial);


  Figure2D fig_bc_image ("box_chains_image", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_bc_image.set_window_properties({650,50},{500,500});
  fig_bc_image.set_axes(axis_limit_image);

  Figure2D fig_intersect ("box_chains_intersections", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_intersect.set_window_properties({1250,50},{500,500});
  fig_intersect.set_axes(axis_limit_image);

  for (size_t i = 0; i < L_BC.size(); ++i)
  {
    const auto& chain = L_BC[i];
    for (const auto& link : chain)
    {
      AnalyticFunction gi ({X},f(link.psi(X)));
      double color_value = (double)i/(double)(L_BC.size()-1);
      fig_bc_init.draw_box(link.psi.eval(link.x), {cmap_edge.color(color_value), cmap_fill.color(color_value)});
      fig_bc_image.draw_box(gi.eval(link.x), {cmap_edge.color(color_value), cmap_fill.color(color_value)});
    }
  }
}
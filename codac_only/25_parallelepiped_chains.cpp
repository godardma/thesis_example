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
        {
          if (!chain.is_neighbor(link))
            cout << "Warning: link is not neighbor of chain" << endl;
          chain.push_back(link);
        }
      else
        break;
    }
    L_BC.push_back(chain);
  }
  return L_BC;
}

bool intersect(const Parallelepiped& p1, const Parallelepiped& p2)
{
  Parallelepiped p_i (p1.c-p2.c, p1.A-p2.A);
  return p_i.box().contains(Vector::constant(p1.c.size(),0));
}

bool find_intersection(const AnalyticFunction<VectorType>& f, vector<ParallelepipedChain>& L_BC)
{
  bool found_intersection = false;
  for (size_t i = 0; i < L_BC.size(); ++i)
  {
    if (i % 10 == 0)
      cout << "Checking chain " << i << " of " << L_BC.size() << endl;
    auto& chain1 = L_BC[i];
    VectorVar X(chain1[0].x.size());
    AnalyticFunction g1 ({X},f(chain1[0].psi(X)));
    for (size_t j = i+1; j < L_BC.size(); ++j)
    {
      auto& chain2 = L_BC[j];
      AnalyticFunction g2 ({X},f(chain2[0].psi(X)));
      for (auto& link1 : chain1)
      {
        for (auto& link2 : chain2)
        {
          auto z1 = g1.parallelepiped_eval(link1.x);
          auto z2 = g2.parallelepiped_eval(link2.x);
          if (intersect(z1,z2))
          {
            auto y1 = link1.psi.eval(link1.x);
            auto y2 = link2.psi.eval(link2.x);
            if (!injectivity_criterion(f, y1|y2))
            {
              found_intersection = true;
              link1.is_intersected = true;
              link2.is_intersected = true;
            }
          }
        }
      }
    }
  }
  return found_intersection;
}


vector<IntervalVector> subdivide(const IntervalVector& domain,
                                        int nx, int ny)
{
    vector<IntervalVector> boxes;

    for (int i = 0; i < nx; i++)
    {
        bool reverse = (i % 2 == 1);

        for (int k = 0; k < ny; k++)
        {
            int j = reverse ? ny - 1 - k : k;

            double x0 = domain[0].lb() +
                        (domain[0].diam() * i) / nx;

            double x1 = domain[0].lb() +
                        (domain[0].diam() * (i+1)) / nx;

            double y0 = domain[1].lb() +
                        (domain[1].diam() * j) / ny;

            double y1 = domain[1].lb() +
                        (domain[1].diam() * (j+1)) / ny;

            boxes.push_back(
                IntervalVector{
                    Interval(x0, x1),
                    Interval(y0, y1)
                }
            );
        }
    }

    return boxes;
}

int main()
{
  // Graphical setup

  ColorMap cmap = ColorMap::rainbow(0.5);


  Figure3D fig_torus("init_torus");
  fig_torus.draw_axes(0.7,{-1.2,-1.2,-1.2});

  Figure3D fig_image("image_torus");
  fig_image.draw_axes(1.,{-2,0,-1});

  // Constructing covers

  IntervalVector V1 ({{-1,1},{-1,1}});
  IntervalVector V2 ({{-1,1},{-1,1}});
  IntervalVector V3 ({{-1,1},{-1,1}});
  IntervalVector V4 ({{-1,1},{-1,1}});
  IntervalVector V5 ({{-1,1},{-1,1}});
  IntervalVector V6 ({{-1,1},{-1,1}});

  auto C1 =subdivide(V1, 20, 20);
  auto C2 =subdivide(V2, 20, 20);
  auto C3 =subdivide(V3, 20, 20);
  auto C4 =subdivide(V4, 20, 20);
  auto C5 =subdivide(V5, 20, 20);
  auto C6 =subdivide(V6, 20, 20);

  cout << "Number of boxes in C1: " << C1.size() << endl;

  // Box Chain decomposition

  double r = 0.8;
  double R = 1.;

  VectorVar X(2);

  AnalyticFunction psi_0 ({X},{1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

  OctaSym id ({1, 2, 3});
  OctaSym s1 ({-2, 1, 3});
  OctaSym s2 ({3, 2, -1});

  AnalyticFunction psi_1 ({X},id(psi_0(X)));
  AnalyticFunction psi_2 ({X},s1(psi_0(X)));
  AnalyticFunction psi_3 ({X},(s1*s1)(psi_0(X)));
  AnalyticFunction psi_4 ({X},(s1.invert())(psi_0(X)));
  AnalyticFunction psi_5 ({X},s2(psi_0(X)));
  AnalyticFunction psi_6 ({X},(s2.invert())(psi_0(X)));

  VectorVar Y (3);
  AnalyticFunction f ({Y},{sqr(Y[0])-sqr(Y[1])+Y[0],2*Y[0]*Y[1]+Y[1],Y[2]});

  AnalyticFunction g1 ({X},f(psi_1(X)));
  AnalyticFunction g2 ({X},f(psi_2(X)));
  AnalyticFunction g3 ({X},f(psi_3(X)));
  AnalyticFunction g4 ({X},f(psi_4(X)));
  AnalyticFunction g5 ({X},f(psi_5(X)));
  AnalyticFunction g6 ({X},f(psi_6(X)));


  for (const auto&box : C1)
  {
    fig_image.draw_parallelepiped(g1.parallelepiped_eval(box), StyleProperties(Color::red(0.5),"atlas"));
    fig_torus.draw_parallelepiped(psi_1.parallelepiped_eval(box), StyleProperties(Color::red(0.5),"atlas"));
  }

  for (const auto&box : C2)
  {
    fig_image.draw_parallelepiped(g2.parallelepiped_eval(box), StyleProperties(Color::green(0.5),"atlas"));
    fig_torus.draw_parallelepiped(psi_2.parallelepiped_eval(box), StyleProperties(Color::green(0.5),"atlas"));
  }

  for (const auto&box : C3)
  {
    fig_image.draw_parallelepiped(g3.parallelepiped_eval(box), StyleProperties(Color::blue(0.5),"atlas"));
    fig_torus.draw_parallelepiped(psi_3.parallelepiped_eval(box), StyleProperties(Color::blue(0.5),"atlas"));
  }

  for (const auto&box : C4)
  {
    fig_image.draw_parallelepiped(g4.parallelepiped_eval(box), StyleProperties(Color::orange(0.5),"atlas"));
    fig_torus.draw_parallelepiped(psi_4.parallelepiped_eval(box), StyleProperties(Color::orange(0.5),"atlas"));
  }

  for (const auto&box : C5)
  {
    fig_image.draw_parallelepiped(g5.parallelepiped_eval(box), StyleProperties(Color::purple(0.5),"atlas"));
    fig_torus.draw_parallelepiped(psi_5.parallelepiped_eval(box), StyleProperties(Color::purple(0.5),"atlas"));
  }

  for (const auto&box : C6)
  {
    fig_image.draw_parallelepiped(g6.parallelepiped_eval(box), StyleProperties(Color::cyan(0.5),"atlas"));
    fig_torus.draw_parallelepiped(psi_6.parallelepiped_eval(box), StyleProperties(Color::cyan(0.5),"atlas"));
  }


  auto L_BC1 = BC_decomposition(f, psi_1, C1);
  auto L_BC2 = BC_decomposition(f, psi_2, C2);
  auto L_BC3 = BC_decomposition(f, psi_3, C3);
  auto L_BC4 = BC_decomposition(f, psi_4, C4);
  auto L_BC5 = BC_decomposition(f, psi_5, C5);
  auto L_BC6 = BC_decomposition(f, psi_6, C6);


  auto L_BC = L_BC1;
  for (const auto& chain : L_BC2)
    L_BC.push_back(chain);
  for (const auto& chain : L_BC3)
    L_BC.push_back(chain);
  for (const auto& chain : L_BC4)
    L_BC.push_back(chain);
  for (const auto& chain : L_BC5)
    L_BC.push_back(chain);
  for (const auto& chain : L_BC6)
    L_BC.push_back(chain);

  cout << "Number of box chains: " << L_BC.size() << endl;

  srand(151);

  for (size_t i = 0; i < L_BC.size(); ++i)
  {
    Color color = Color::random(0.5);
    const auto& chain = L_BC[i];
    for (const auto& link : chain)
    {
      AnalyticFunction gi ({X},f(link.psi(X)));
      fig_image.draw_parallelepiped(gi.parallelepiped_eval(link.x), StyleProperties(color,"box_chains"));
      fig_torus.draw_parallelepiped(link.psi.parallelepiped_eval(link.x), StyleProperties(color,"box_chains"));
    }
  }

  // Detecting self-intersections
  cout << "looking for intersections..." << endl;

  bool intersection_found = find_intersection(f, L_BC);
  cout << "intersection found ? "<< (intersection_found?"true":"false") << endl;

  vector<ParallelepipedLink> unintersected_links;
  vector<ParallelepipedLink> intersected_links;
  vector<ParallelepipedLink> L_BL;

  for (const auto& chain : L_BC)
  {
    for (const auto& link : chain)
    {
      L_BL.push_back(link);
      AnalyticFunction gi ({X},f(link.psi(X)));
      if (link.is_intersected)
        {
          fig_image.draw_parallelepiped(gi.parallelepiped_eval(link.x), StyleProperties(Color::red(0.5),"intersected"));
          fig_torus.draw_parallelepiped(link.psi.parallelepiped_eval(link.x), StyleProperties(Color::red(0.5),"intersected"));
        }
      else
      {
        fig_image.draw_parallelepiped(gi.parallelepiped_eval(link.x), StyleProperties(Color::gray(0.5),"intersected"));
        fig_torus.draw_parallelepiped(link.psi.parallelepiped_eval(link.x), StyleProperties(Color::gray(0.5),"intersected"));
      }
    }
  }
}
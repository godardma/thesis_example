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

int main()
{
  IntervalVector V1 ({{-1,1}});
  IntervalVector V2 ({{-1,1}});

  vector<IntervalVector> C1;
  vector<IntervalVector> C2;

  double epsilon = 0.1;

  split(V1, epsilon, C1);
  split(V2, epsilon, C2);

  VectorVar X(1);
  AnalyticFunction psi_1 ({X},{cos(X[0]*PI/2),sin(X[0]*PI/2)});
  AnalyticFunction psi_2 ({X},{cos(PI+X[0]*PI/2),sin(PI+X[0]*PI/2)});

  VectorVar Y(2);
  AnalyticFunction f ({Y},{sqr(Y[0])-sqr(Y[1])+Y[0],2*Y[0]*Y[1]+Y[1]});

  AnalyticFunction g1 ({X},f(psi_1(X)));
  AnalyticFunction g2 ({X},f(psi_2(X)));

  cout<<injectivity_criterion(f, chain1, link1)<<endl;;

  for (auto& c1 : C1)
    DefaultFigure::draw_box(g1.eval(c1));

  for (auto& c2 : C2)
    DefaultFigure::draw_box(g2.eval(c2));
}
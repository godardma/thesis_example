#include <codac>
#include <vector>

using namespace codac2;
using namespace std;

class BoxLink
{
  public:
    AnalyticFunction<VectorType> psi;
    IntervalVector x;
    bool is_intersected = false;
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

bool injectivity_criterion(const AnalyticFunction<VectorType>& f, const IntervalVector& x)
{
  IntervalMatrix Jf = f.diff(x);
  IntvFullPivLU lu(Jf);
  return !lu.determinant().contains(0);
}

bool injectivity_criterion (const AnalyticFunction<VectorType>& f, const BoxChain& chain, const BoxLink& link)
{
  return injectivity_criterion(f,link.psi.eval(chain.bounding_box() | link.x));
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

bool find_intersection(const AnalyticFunction<VectorType>& f, vector<BoxChain>& L_BC)
{
  bool found_intersection = false;
  for (size_t i = 0; i < L_BC.size(); ++i)
  {
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
          auto z1 = g1.eval(link1.x);
          auto z2 = g2.eval(link2.x);
          if (z1.intersects(z2))
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

vector<BoxChain> extract_new_chains(const AnalyticFunction<VectorType>& f,const vector<BoxLink>& L_U)
{
  VectorVar X (L_U[0].x.size());
  vector<BoxChain> L_BC;
  vector<const BoxLink*> L_w;
  for (const auto& link : L_U)
      L_w.push_back(&link);
  BoxChain L_B;
  bool was_added = false;
  int nb_added = 0;
  while (!L_w.empty())
  {
    if (!was_added)
    {
      L_B.push_back(*L_w[0]);
      L_w.erase(L_w.begin());
    } 
    was_added = false;
    for (size_t i = 0; i < L_w.size(); ++i)
    {
      AnalyticFunction g_w ({X},f(L_w[i]->psi(X)));
      for (const auto& link : L_B)
      {
        AnalyticFunction g_link ({X},f(link.psi(X)));
        if (g_w.eval(L_w[i]->x).intersects(g_link.eval(link.x))) // we can add it to the current box chain
        {
          L_B.push_back(*L_w[i]);
          L_w.erase(L_w.begin()+i);
          was_added = true;
          nb_added++;
          break;
        }
      }
      if (was_added)
        break;      
    }
    if (!was_added)
    {
      L_BC.push_back(L_B);
      L_B.clear();
    }
  }
  if (!L_B.empty())
    L_BC.push_back(L_B);
  return L_BC;
}

PavingOut pave_from_boundary(const AnalyticFunction<VectorType>& f, const vector<BoxLink>& L_BL, const IntervalVector& domain, double epsilon)
{
  CtcUnion ctc_union(L_BL[0].x.size());
  VectorVar X(L_BL[0].x.size());

  for (const auto& link : L_BL)
  {
    AnalyticFunction g ({X},f(link.psi(X)));
    IntervalVector z = g.eval(link.x);
    ctc_union |= CtcWrapper(z);
  }
  return pave(domain, ctc_union, epsilon);
}

vector<PavingOut::ConnectedSubset_> get_inner_areas(const AnalyticFunction<VectorType>& f, const PavingOut& p, const vector<BoxLink>& L_BL, double delta = 0.1)
{
  vector<PavingOut::ConnectedSubset_> inner_areas;
  auto v_cs = p.connected_subsets(PavingOut::outer_complem);

  for (const auto& cs : v_cs)
  {
    bool is_inner = false;
    for (const auto& link : L_BL)
    {
      VectorVar X(link.x.size());
      AnalyticFunction g ({X},f(link.psi(X)));
      IntervalVector z = g.eval(link.x);
      IntervalMatrix Jg = g.diff(link.x);
      Vector tangent = Jg.mid().col(0);
      Vector normal ({-tangent[1], tangent[0]});
      Vector point_in = z.mid() + delta*normal;
      for (const auto& box : cs.boxes())
      {
        if (box.contains(point_in))
        {
          inner_areas.push_back(cs);
          is_inner = true;
          break;
        }
      }
      if (is_inner)
        break;
    }
  }
  return inner_areas;
}

vector<BoxChain> remove_fake(const AnalyticFunction<VectorType>& f, const vector<BoxChain>& L_BC, const vector<PavingOut::ConnectedSubset_>& inner_areas, double delta = 0.1)
{
  vector<BoxChain> L_BC_no_fake;
  for (const auto& chain : L_BC)
  {
    bool is_fake = false;
    for (const auto& link : chain)
    {
      VectorVar X(link.x.size());
      AnalyticFunction g ({X},f(link.psi(X)));
      IntervalVector z = g.eval(link.x);
      IntervalMatrix Jg = g.diff(link.x);
      Vector tangent = Jg.mid().col(0);
      Vector normal ({-tangent[1], tangent[0]});
      Vector point_out = z.mid() - delta*normal;
      for (const auto& cs : inner_areas)
      {
        for (const auto& box : cs.boxes())
        {
          if (box.contains(point_out))
          {
            is_fake = true;
            break;
          }
        }
        if (is_fake)
          break;
      }
    }
    if (!is_fake)
      L_BC_no_fake.push_back(chain);

  }
  return L_BC_no_fake;
}

int main()
{
  // Graphical setup

  ColorMap cmap_edge = ColorMap::rainbow();
  ColorMap cmap_fill = ColorMap::rainbow(0.5);

  IntervalVector axis_limit_initial ({{-1.5,1.5},{-1.5,1.5}});
  IntervalVector axis_limit_image ({{-2.25,2.75},{-2.5,2.5}});

  Figure2D fig_box_chains_new("box_chains_new", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_box_chains_new.set_window_properties({350,50},{500,500});
  fig_box_chains_new.set_axes(axis_limit_image);

  Figure2D fig_fake_boundary_removed("fake_boundary_removed", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_fake_boundary_removed.set_window_properties({950,50},{500,500});
  fig_fake_boundary_removed.set_axes(axis_limit_image);

  Figure2D fig_paving_new_in_out("paving_new_in_out", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_paving_new_in_out.set_window_properties({350,500},{500,500});
  fig_paving_new_in_out.set_axes(axis_limit_image);

  Figure2D fig_bc_init("box_chains_init", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_bc_init.set_window_properties({50,50},{500,500});
  fig_bc_init.set_axes(axis_limit_initial);

  Figure2D fig_bc_image ("box_chains_image", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_bc_image.set_window_properties({650,50},{500,500});
  fig_bc_image.set_axes(axis_limit_image);
  
  Figure2D fig_intersect ("box_chains_intersections", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_intersect.set_window_properties({1250,50},{500,500});
  fig_intersect.set_axes(axis_limit_image);

  Figure2D fig_paving_out("paving_out", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_paving_out.set_window_properties({50,650},{500,500});
  fig_paving_out.set_axes(axis_limit_image);

  Figure2D fig_paving_out_vcs("paving_out_vcs", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_paving_out_vcs.set_window_properties({650,650},{500,500});
  fig_paving_out_vcs.set_axes(axis_limit_image);

  Figure2D fig_paving_in_out("paving_in_out", GraphicOutput::VIBES | GraphicOutput::IPE);
  fig_paving_in_out.set_window_properties({1250,650},{500,500});
  fig_paving_in_out.set_axes(axis_limit_image);

  // Constructing covers

  IntervalVector V1 ({{-1,1}});
  IntervalVector V2 ({{-1,1}});

  vector<IntervalVector> C1;
  vector<IntervalVector> C2;

  double epsilon = 0.1;

  split(V1, epsilon, C1);
  split(V2, epsilon, C2);

  cout << "Number of boxes in cover: " << C1.size() << endl;

  // Box Chain decomposition

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

  // Detecting self-intersections

  bool intersection_found = find_intersection(f, L_BC);
  cout << "intersection found ? "<< (intersection_found?"true":"false") << endl;

  vector<BoxLink> unintersected_links;
  vector<BoxLink> intersected_links;
  vector<BoxLink> L_BL;

  for (const auto& chain : L_BC)
  {
    for (const auto& link : chain)
    {
      L_BL.push_back(link);
      AnalyticFunction gi ({X},f(link.psi(X)));
      if (link.is_intersected)
      {
        intersected_links.push_back(link);
        fig_box_chains_new.draw_box(gi.eval(link.x), {Color::red(), Color::red(0.5)});
        fig_intersect.draw_box(gi.eval(link.x), {Color::red(), Color::red(0.5)});
        fig_fake_boundary_removed.draw_box(gi.eval(link.x), {Color::red(), Color::red(0.5)});

      }
      else
      {
        unintersected_links.push_back(link);
        fig_intersect.draw_box(gi.eval(link.x), {Color::black(), Color::black(0.5)});
      }
    }
  }

  // Extracting new box chains from unintersected links

  auto L_BC_new = extract_new_chains(f, unintersected_links);

  for (size_t i = 0; i < L_BC_new.size(); ++i)
  {
    const auto& chain = L_BC_new[i];
    for (const auto& link : chain)
    {
      AnalyticFunction gi ({X},f(link.psi(X)));
      double color_value = (double)(&chain-&L_BC_new[0])/(double)(L_BC_new.size()-1);
      if (i==0)
        fig_box_chains_new.draw_box(gi.eval(link.x), {Color::blue(), Color::blue(0.5)});
      else if (i==1)
        fig_box_chains_new.draw_box(gi.eval(link.x), {Color::green(), Color::green(0.5)});
    }
  }

  // Initial paving

  auto p = pave_from_boundary(f, L_BL, axis_limit_image, epsilon);
  fig_paving_out.draw_paving(p);

  // Paving connected subsets

  fig_paving_out_vcs.draw_paving(p);
  auto v_cs = p.connected_subsets(PavingOut::outer_complem);
  int i = 0;
  for (const auto& cs : v_cs)
  {
    if (i==0)
      fig_paving_out_vcs.draw_subpaving(cs,{Color::black(),Color::red(0.8)});
    else if (i==1)
      fig_paving_out_vcs.draw_subpaving(cs,{Color::black(),Color::green(0.8)});
    else if (i==2)
      fig_paving_out_vcs.draw_subpaving(cs,{Color::black(),Color::blue(0.8)});
    ++i;
  }

  // Determining inner areas

  auto inner_areas = get_inner_areas(f, p, unintersected_links);

  fig_paving_in_out.draw_paving(p);
  for (const auto& cs : inner_areas)
    fig_paving_in_out.draw_subpaving(cs,StyleProperties::inside());
  
  // Removing fake box chains

  auto L_BC_no_fake = remove_fake(f, L_BC_new, inner_areas);
  vector<BoxLink> L_no_fake;

  for (const auto& chain : L_BC_no_fake)
  {
    for (const auto& link : chain)
    {
      AnalyticFunction gi ({X},f(link.psi(X)));
      fig_fake_boundary_removed.draw_box(gi.eval(link.x), {Color::black(), Color::black(0.5)});
      L_no_fake.push_back(link);
    }
  }

  // Adding intersected box as part of the real boundary

  for (const auto& link : intersected_links)
  {
    AnalyticFunction gi ({X},f(link.psi(X)));
    L_no_fake.push_back(link);
  }

  // New paving from the new boundary

  auto p_no_fake = pave_from_boundary(f, L_no_fake, axis_limit_image, epsilon);
  auto inner_areas_no_fake = get_inner_areas(f, p_no_fake, L_no_fake);
  fig_paving_new_in_out.draw_paving(p_no_fake);
  for (const auto& cs : inner_areas_no_fake)
    fig_paving_new_in_out.draw_subpaving(cs,StyleProperties::inside());

}
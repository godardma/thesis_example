#include <codac>
#include <codac-capd.h>

#include "wrappers_utils.h"

using namespace std;
using namespace codac2;

class ZonotopeLink
{
public:
    AnalyticFunction<VectorType> psi;
    IntervalVector x;
    IntervalMatrix Jf;
    Zonotope z;
    Zonotope q;
    bool is_intersected = false;
    ZonotopeLink(const AnalyticFunction<VectorType> &psi, const IntervalVector &x, Zonotope &z, Zonotope &q, const IntervalMatrix &Jf) : psi(psi), x(x), z(z), q(q), Jf(Jf) {}

    bool is_neighbor(const ZonotopeLink &other) const
    {
        return x.intersects(other.x);
    }
};

class ZonotopeChain : public vector<ZonotopeLink>
{
public:
    bool is_neighbor(const ZonotopeLink &link)
    {
        for (const auto &l : *this)
            if (l.is_neighbor(link))
                return true;
        return false;
    }

    IntervalVector bounding_box() const
    {
        IntervalVector box = (*this)[0].x;
        for (const auto &link : *this)
            box |= link.x;
        return box;
    }
};

bool injectivity_criterion(const ZonotopeChain &chain, const ZonotopeLink &link)
{
    auto Jphi = link.Jf;
    for (const auto &chain_link : chain)
        Jphi |= chain_link.Jf;

    IntvFullPivLU lu(Jphi);
    return !lu.determinant().contains(0);
}

bool injectivity_criterion_recheck(const capd::IMap &IMap, double tf, const ZonotopeLink &link1, const ZonotopeLink &link2)
{
    IntervalVector y1 = link1.psi.eval(link1.x);
    IntervalVector y2 = link2.psi.eval(link2.x);

    capd::IMap g(IMap);
    capd::IOdeSolver solver(g, 30);

    solver.setAbsoluteTolerance(1e-20);
    solver.setRelativeTolerance(1e-20);

    capd::ITimeMap timeMap(solver);
    capd::ITimeMap timeMap_punc(solver);

    capd::interval initialTime(0.);
    capd::interval finalTime(tf);

    IntervalVector Y = y1|y2;

    try
    {
        capd::IMatrix monodromyMatrix(4,4);
        capd::ITimeMap::SolutionCurve solution(initialTime);
        capd::IVector c = to_capd(Y);

        capd::C1Rect2Set s(c);
        timeMap(finalTime, s, solution);
        capd::IVector result = timeMap(finalTime, s, monodromyMatrix);
        IntervalMatrix JJf = to_codac(monodromyMatrix);

        IntervalMatrix A = JJf.block(0, 0, 2, 2);
        IntervalMatrix B = JJf.block(0, 2, 2, 2);

        double rad = 0.4;
        IntervalMatrix Jphi = A + rad * B;

        IntvFullPivLU lu(Jphi);

        return !lu.determinant().contains(0);
    }
    catch(...)
    {
        return false;
    }

}


ZonotopeLink PEIBOS3D_step(const capd::IMap &IMap, double tf, IntervalVector X, const AnalyticFunction<VectorType> &psi_0, const OctaSym &symmetry, double epsilon, const Vector &offset)
{
    int m = psi_0.input_size();
    int n = psi_0.output_size();

    assert(offset.size() == n);
    assert(m < n);
    assert(symmetries.size() > 0 && (int)symmetries[0].size() == n);

    clock_t t_start = clock();

    VectorVar X_var(m);
    AnalyticFunction psi({X_var}, symmetry(psi_0(X_var)));

    vector<Parallelepiped> output;

    // CAPD solver setup
    capd::IMap g(IMap);
    capd::IOdeSolver solver(g, 30);

    solver.setAbsoluteTolerance(1e-20);
    solver.setRelativeTolerance(1e-20);

    capd::ITimeMap timeMap(solver);
    capd::ITimeMap timeMap_punc(solver);

    capd::interval initialTime(0.);
    capd::interval finalTime(tf);

    // To get the flow function and its Jacobian (monodromy matrix) for [x]
    IntervalVector Y = psi.eval(X) + offset;

    capd::IMatrix monodromyMatrix(n, n);
    capd::ITimeMap::SolutionCurve solution(initialTime);
    capd::IVector c = to_capd(Y);

    capd::C1Rect2Set s(c);
    timeMap(finalTime, s, solution);
    capd::IVector result = timeMap(finalTime, s, monodromyMatrix);
    IntervalMatrix JJf = to_codac(monodromyMatrix);

    // To get the flow function and its Jacobian (monodromy matrix) for x_hat
    Vector xc = X.mid();
    IntervalVector yc = (psi.eval(xc) + offset);

    capd::IMatrix monodromyMatrix_punc(n, n);
    capd::ITimeMap::SolutionCurve solution_punct(initialTime);
    capd::IVector c_punct = to_capd(IntervalVector(yc));

    capd::C1Rect2Set s_punct(c_punct);
    timeMap_punc(finalTime, s_punct, solution_punct);
    capd::IVector result_punct = timeMap_punc(finalTime, s_punct, monodromyMatrix_punc);
    IntervalMatrix JJf_point = to_codac(monodromyMatrix_punc);

    // Center of the parallelepiped
    Vector z = Vector(to_codac(result).mid());

    auto p = parallelepiped_inclusion(to_codac(result_punct), JJf, JJf_point.mid(), psi_0, symmetry, X);
    IntervalMatrix A = JJf.block(0, 0, 2, 2);
    IntervalMatrix B = JJf.block(0, 2, 2, 2);

    double rad = 0.4;
    IntervalMatrix Jphi = A + rad * B;

    Zonotope zonotope = p.proj({0, 1});
    Zonotope q = p.proj({2, 3});
    ZonotopeLink zlink(psi, X, zonotope,q, Jphi);
    return zlink;
}

template <typename T>
vector<ZonotopeChain> PEIBOS_custom(const capd::IMap &i_map_wrapper, double tf, AnalyticFunction<T> &psi_0, const vector<OctaSym> &symmetries, double epsilon, Vector offset)
{
    vector<IntervalVector> cover;
    for (double t1 = -1; t1 < 1; t1 += epsilon)
        cover.push_back(IntervalVector({{t1, t1 + epsilon}}));
    cout << cover.size() << endl;
    vector<ZonotopeChain> v_z_chain;
    for (int i = 0; i < symmetries.size(); i++)
    {
        auto symmetry = symmetries[i];

        auto L_W = cover;
        IntervalVector x = L_W.back();
        L_W.pop_back();
        while (!L_W.empty())
        {
            auto zlink = PEIBOS3D_step(i_map_wrapper, tf, x, psi_0, symmetry, epsilon, offset);
            ZonotopeChain z_chain;
            z_chain.push_back(zlink);
            while (!L_W.empty())
            {
                x = L_W.back();
                L_W.pop_back();
                auto zlink = PEIBOS3D_step(i_map_wrapper, tf, x, psi_0, symmetry, epsilon, offset);
                if (injectivity_criterion(z_chain,zlink) && z_chain.is_neighbor(zlink))
                        z_chain.push_back(zlink);
                else
                    break;
            }
            v_z_chain.push_back(z_chain);
        }
    }
    return v_z_chain;
}

bool find_intersection(const capd::IMap &i_map_wrapper, double tf, vector<ZonotopeChain>& L_BC)
{
  bool found_intersection = false;
  for (size_t i = 0; i < L_BC.size(); ++i)
  {
    auto& chain1 = L_BC[i];
    for (size_t j = i+1; j < L_BC.size(); ++j)
    {
      auto& chain2 = L_BC[j];
      for (auto& link1 : chain1)
      {
        for (auto& link2 : chain2)
        {
          if ((link1.z.box()).intersects(link2.z.box()))
          {
            IntervalVector y1 = link1.psi.eval(link1.x);
            IntervalVector y2 = link2.psi.eval(link2.x);
            if ((y1&y2).is_empty())
            {
                if (!injectivity_criterion_recheck(i_map_wrapper, tf, link1,link2))
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
  }
  return found_intersection;
}

int main()
{
    set_nb_threads(max_threads());

    auto wrappers = readIMapWrappers("../wrappers/cshape_wrappers.txt");

    double rad = 0.4;
    double x1_init = 0., x2_init = 0.0;

    capd::IMap vectorField_wrap(std::get<2>(wrappers));

    double tf_discrete = 8.;
    double dt_discrete = 8.;

    IntervalVector Y0 ({{3, 11},{-1.5, 1.5}});

    Figure2D output_discrete("cshape_init_chains", GraphicOutput::VIBES | GraphicOutput::IPE);
    output_discrete.set_axes(axis(0, {-0.5, 0.5}), axis(1, {-0.5, 0.5}));
    output_discrete.set_window_properties({50, 200}, {800, 800});

    Figure2D output_final("cshape_chains", GraphicOutput::VIBES | GraphicOutput::IPE);
    output_final.set_axes(Y0);
    output_final.set_window_properties({50, 200}, {1600, 800});

    Figure2D output_intersections("cshape_chains_intersections", GraphicOutput::VIBES | GraphicOutput::IPE);
    output_intersections.set_axes(Y0);
    output_intersections.set_window_properties({50, 200}, {1600, 800});

    Figure2D output_pave_in_out("cshape_chains_inout", GraphicOutput::IPE);
    output_pave_in_out.set_axes(Y0);
    output_pave_in_out.set_window_properties({50, 200}, {1600, 800});

    Figure2D output_cleaned("cshape_cleaned", GraphicOutput::VIBES | GraphicOutput::IPE);
    output_cleaned.set_axes(Y0);
    output_cleaned.set_window_properties({50, 200}, {1600, 800});

    // Reachability analysis

    cout << "starting reachability" << endl;
    VectorVar X_reach(1);
    // AnalyticFunction psi0_reach({X_reach}, {rad * cos(X_reach[0] * PI / 2.5), rad * sin(X_reach[0] * PI / 2.5), cos(X_reach[0] * PI / 2.5), sin(X_reach[0] * PI / 2.5)});
    AnalyticFunction psi0_reach({X_reach}, {rad * cos(X_reach[0] * PI / 4.), rad * sin(X_reach[0] * PI / 4.), cos(X_reach[0] * PI / 4.), sin(X_reach[0] * PI / 4.)});

    OctaSym id_reach({1, 2, 3, 4});
    OctaSym s({-2, 1, -4, 3});

    vector<OctaSym> Sigma({id_reach, s, s * s, s.invert()});
    // vector<double> resolutions ({0.01,0.001,0.01,0.001});
    // vector<double> resolutions ({0.01,0.00025,0.002,0.00025});
    vector<double> resolutions ({0.005,0.00025,0.002,0.00025});



    ColorMap cmap_edge = ColorMap::rainbow();
    ColorMap cmap_fill = ColorMap::rainbow(0.5);

    vector<ZonotopeChain> v_z_chain;

    for (int i =0; i<Sigma.size();i++)
    {
        vector<ZonotopeChain> v_z_chain_i = PEIBOS_custom(vectorField_wrap, tf_discrete, psi0_reach, {Sigma[i]}, resolutions[i], {x1_init, x2_init, 0, 0});
        double max_i = ((double) v_z_chain_i.size()) -1.0;
        double index = 0;
        for (const auto &z_chain : v_z_chain_i)
        {
            v_z_chain.push_back(z_chain);
            for (const auto &z_link : z_chain)
            {
                output_final.draw_zonotope(z_link.z, StyleProperties({cmap_edge.color(index/max_i), cmap_fill.color(index/max_i)}, "reachable_final", "z:1"));
                auto p = z_link.psi.parallelepiped_eval(z_link.x);
                output_discrete.draw_zonotope(p.proj({0,1}), StyleProperties({cmap_edge.color(index/max_i), cmap_fill.color(index/max_i)}, "reachable_final", "z:1"));

            }
            index+=1.0;
        }
    }

    cout << v_z_chain.size() << endl;
    find_intersection(vectorField_wrap, tf_discrete,v_z_chain);

    
    CtcUnion ctc_union(2);

    vector<ZonotopeLink> unintersected_links;
    vector<ZonotopeLink> kept_links;
    
    for (const auto& chain : v_z_chain)
    {
        for (const auto & link:chain)
        {

            if (link.is_intersected)
            {
                kept_links.push_back(link);
                output_intersections.draw_zonotope(link.z, StyleProperties({Color::red(),  Color::red(0.5)},"z:1"));           
            }
            else
            {
                output_intersections.draw_zonotope(link.z, StyleProperties({Color::gray(),  Color::gray(0.5)})); 
                unintersected_links.push_back(link);
          
            }
            ctc_union |= CtcWrapper(link.z.box());
        }
    }

    int dim = 2;
    SepCtcPair sep_boundary(CtcIdentity(dim), ctc_union);

    PavingInOut p = pave(Y0,sep_boundary,0.01);

    auto v_cs = p.connected_subsets(PavingInOut::outer_complem);
    std::list<std::shared_ptr<PavingInOut::ConnectedSubset_>> cs_to_color;

    double eps_color = 0.01;

    for (const auto& cs: v_cs)
    {
        bool to_color = false;
        for (const auto box : cs.boxes())
        {
            if (to_color)
            {
                cs_to_color.push_back(std::make_shared<PavingInOut::ConnectedSubset_>(cs));
                break;
            }
            for (const auto link : unintersected_links)
            {
                Vector behind = link.z.c-eps_color*link.q.c;
                DefaultFigure::draw_point(behind);
                if ((box).contains(behind))
                {
                    DefaultFigure::draw_point(behind,StyleProperties({Color::blue(),Color::blue()},"z:1"));  
                    to_color = true;
                    break;
                }
            }
        }
    }

    auto visitor = [&](std::shared_ptr<PavingNode<PavingInOut>> n)
    {
        std::list<IntervalVector> boxes_to_color;
        for (const auto& cs : cs_to_color)
            for (const auto& box : cs->boxes())
            boxes_to_color.push_back(box);

        IntervalVector h = n->hull();
        IntervalVector u = n->unknown();

        for (const auto& bi : boxes_to_color)
        {
            IntervalVector bih = bi & h;
            if (!bih.is_empty())
            {
            if (std::get<0>(n->boxes()).is_empty())
            { 
                n->set_boxes(std::make_tuple(std::get<1>(n->boxes()),std::get<0>(n->boxes())));
                return false;
            }

            std::list<IntervalVector> lbi;

            for (const auto& bj : boxes_to_color)
            {
                IntervalVector bjh = bj & h;
                if (!bjh.is_empty())
                lbi.push_back(bjh);
            }

            IntervalVector prev_x_in = IntervalVector::empty(dim);

            // there used to be a copy x_in = h there, but didn't behave as exepected

            while (h != prev_x_in)
            {
                prev_x_in = IntervalVector(h);

                for (const auto& li : lbi)
                {
                std::list<IntervalVector> d = h.diff(li);

                if (!d.empty())
                {
                    
                    IntervalVector hi = IntervalVector::empty(dim);
                    for (const auto& di : d)
                    hi |= (h & di);
                    h &= hi;

                }
                }
            }
            
            n->set_boxes(std::make_tuple(hull(lbi)|u, h));
            return true;
            }
        }
        return true;
    };

    // Use the visitor to color the inner regions
    p.tree()->visit(visitor);


    DefaultFigure::set_axes(Y0);

    output_pave_in_out.draw_paving(p);


    auto boxes_in = p.boxes(PavingInOut::inner);

    for (const auto link : unintersected_links)
    {
        bool to_keep = true;
        Vector front = link.z.c+eps_color*link.q.c;
        DefaultFigure::draw_point(front);
        for (const auto& box : boxes_in)
        {
            if ((box).contains(front))
            {
                DefaultFigure::draw_point(front,StyleProperties({Color::green(),Color::green()},"z:1"));
                to_keep = false;
                break;
            }
        }
        if (to_keep)
            kept_links.push_back(link);
    }

    for (const auto link : kept_links)
        output_cleaned.draw_zonotope(link.z, StyleProperties::boundary());

}
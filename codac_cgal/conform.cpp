// Author: Maël GODARD
#include "triangulate_tools.h"

int main()
{  
    VectorVar x(1);
    AnalyticFunction f ({x},{cos(x[0]),sin(x[0])});
    
    vector<Parallelepiped> L_p;
    int nb_parallelepiped = 15;
    double delta_theta = 2.0*PI/(double)nb_parallelepiped;

    CtcUnion ctc_union (2);

    for (int i = 0; i < nb_parallelepiped; ++i)
    {
        IntervalVector X ({Interval(i*delta_theta, (i+1)*delta_theta)});
        Parallelepiped p = f.parallelepiped_eval(X);
        L_p.push_back(p);
        ctc_union |= CtcParallelepiped(p);
    }

    IntervalVector x0({{-1.5,1.5},{-1.5,1.5}});
    auto p = pave(x0, ctc_union, 0.15);

    Figure2D fig1 ("Paving_of_circle_parallelelepiped", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig1.set_window_properties({600,50},{500,500});
    fig1.set_axes(x0);

    Figure2D fig2 ("Paving_of_circle_parallelelepiped_inout", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig2.set_window_properties({1150,50},{500,500});
    fig2.set_axes(x0);

    fig1.draw_paving(p);

    for (const auto& par : L_p)
        fig1.draw_parallelepiped(par, StyleProperties(Color::red(),"w:0.008","--"));

    auto v_cs = p.connected_subsets(PavingOut::outer_complem);

    int i = 0;
    for (auto& cs : v_cs)
    {
        if (i==0)
            fig2.draw_subpaving(cs, StyleProperties::outside());
        else
            fig2.draw_subpaving(cs, StyleProperties::inside());
        i++;
    }

    for (const auto& box : p.boxes(PavingOut::outer))
        fig2.draw_box(box, StyleProperties::boundary());

    triangulate(L_p, x0, "triangulation_of_circle_parallelepiped");
}

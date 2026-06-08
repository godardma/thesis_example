#include <codac>
#include <codac-capd.h>

#include "wrappers_utils.h"

using namespace std;
using namespace codac2;


int main()
{
    set_nb_threads(max_threads());

    auto wrappers = readIMapWrappers("../wrappers/cshape_wrappers.txt");

    double rad = 0.4;
    double x1_init = 0., x2_init = 0.0;

    capd::IMap vectorField_wrap(std::get<2>(wrappers));

    double tf_discrete = 8.;
    double dt_discrete = 4.;

    VectorVar X_2d(1);
    AnalyticFunction psi0_2d({X_2d}, {rad * cos(X_2d[0] * PI / 4.), rad * sin(X_2d[0] * PI / 4.)});

    OctaSym id_2d({1, 2});
    OctaSym s_2d({-2, 1});

    Figure2D output_discrete("cshape", GraphicOutput::VIBES | GraphicOutput::IPE);
    output_discrete.set_axes(axis(0, {-1, 11}), axis(1, {-1.5, 1.5}));
    output_discrete.set_window_properties({50, 200}, {1600, 800});

    Figure2D output_final("cshape_final", GraphicOutput::VIBES | GraphicOutput::IPE);
    output_final.set_axes(axis(0, {3, 11}), axis(1, {-1.5, 1.5}));
    output_final.set_window_properties({50, 200}, {800, 800});

    // Undisturbed case

    capd::IMap vectorField_wrap_2d(std::get<0>(wrappers));

    auto peibos_output_base = PEIBOS(vectorField_wrap_2d, tf_discrete, dt_discrete, psi0_2d, {id_2d, s_2d, s_2d * s_2d, s_2d.invert()}, 0.1, {x1_init, x2_init}, true);


    auto m_v_par_base = reach_set(peibos_output_base);
    for (auto &[t, v_par] : m_v_par_base)
        for (const auto &p : v_par)
            {
                output_discrete.draw_parallelepiped(p, StyleProperties({Color::pink(), Color::pink(0.5)},"basic_at_"+to_string(t),"z:1"));
                if (t == tf_discrete)
                    {
                        output_final.draw_parallelepiped(p, StyleProperties({Color::pink(), Color::pink(0.5)},"basic_final","z:1"));
                    }
            }

    // Disturbed case

    srand(158);

    // generate 50 perturbed vector fields
    for (int i=0; i<50; i++)
    {
        double r = ((double)rand()/(double)RAND_MAX)*0.1;
        double theta = ((double)rand()/(double)RAND_MAX)*2.*PI;
        double dx = r*cos(theta);
        double dy = r*sin(theta);
        auto f_disturbed = std::get<1>(wrappers);
        // replace "w1" and "w2" in f_disturbed with dx and dy
        size_t pos_w1 = f_disturbed.find("w1");
        if (pos_w1 != std::string::npos)
            f_disturbed.replace(pos_w1, 2, std::to_string(dx));
        size_t pos_w2 = f_disturbed.find("w2");
        if (pos_w2 != std::string::npos)
            f_disturbed.replace(pos_w2, 2, std::to_string(dy));

        auto peibos_output_disturbed = PEIBOS(f_disturbed, tf_discrete, dt_discrete, psi0_2d, {id_2d, s_2d, s_2d * s_2d, s_2d.invert()}, 0.05, {x1_init, x2_init});
        
        auto m_v_par_disturb = reach_set(peibos_output_disturbed);
        for (auto &[t, v_par] : m_v_par_disturb)
        {
            for (const auto &p : v_par)
            {
                output_discrete.draw_parallelepiped(p, StyleProperties({Color::gray(), Color::gray(0.5)},"perturbed_at_"+to_string(t)));
                if (t == tf_discrete)
                    {
                        output_final.draw_parallelepiped(p, StyleProperties({Color::gray(), Color::gray(0.5)},"perturbed_final"));
                    }
            }
        }
            
    }

    // Reachability analysis

    cout << "starting reachability" << endl;
    VectorVar X_reach(1);
    AnalyticFunction psi0_reach({X_reach}, {rad * cos(X_reach[0] * PI / 4.), rad * sin(X_reach[0] * PI / 4.), cos(X_reach[0] * PI / 4.), sin(X_reach[0] * PI / 4.)});

    OctaSym id_reach({1, 2, 3, 4});
    OctaSym s({-2, 1, -4, 3});
    

    vector<OctaSym> Sigma ({id_reach,s,s*s,s.invert()});

    for (int i=0; i< Sigma.size(); i++)
    {
        auto sigma = Sigma[i];

        auto peibos_output_reach = PEIBOS(vectorField_wrap,tf_discrete,dt_discrete,psi0_reach, {sigma} , 0.0005, {x1_init,x2_init,0,0}, true);

        auto m_v_par_reach = reach_set(peibos_output_reach);

        for (auto& [t, v_par] : m_v_par_reach)
            for (const auto& p : v_par)
            {
                output_discrete.draw_parallelepiped(p, StyleProperties({Color::red(), Color::red(0.5)},"reachable_at_"+to_string(t),"z:1"));
                if (t == tf_discrete)
                    {
                        output_final.draw_parallelepiped(p, StyleProperties({Color::red(), Color::red(0.5)},"reachable_final","z:1"));
                    }
            }
    }
}
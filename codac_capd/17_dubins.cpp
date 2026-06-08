#include <codac>
#include <codac-capd.h>

#include "wrappers_utils.h"

using namespace std;
using namespace codac2;

int main()
{
    SampledTraj<Vector> state_traj;
    std::ifstream in("../datasets/heading_following_west_1.cdc", std::ios::binary);
    deserialize(in, state_traj);

    double tf_discrete = 1.83082;
    double dt_discrete = 0.5;

    Figure3D output("dubins_traj");
    output.draw_axes(1,{-4,0,0});

    Interval x_bounds (-8.3,-7.1);
    Interval y_bounds (-0.25,0.65);
    Interval z_bounds (0.15,0.45);

    Interval traj_x_bounds (-8.5,-4.5);
    Interval traj_y_bounds (-1,1);
    Interval traj_z_bounds (-0.2,0.5);


    Figure2D fig_traj_xy("traj_xy", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig_traj_xy.set_window_properties({50, 50}, {500, 500});
    fig_traj_xy.set_axes(IntervalVector({traj_x_bounds,traj_y_bounds}));

    Figure2D fig_traj_yz("traj_yz", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig_traj_yz.set_window_properties({600, 50}, {500, 500});
    fig_traj_yz.set_axes(IntervalVector({traj_y_bounds,traj_z_bounds}));

    Figure2D fig_traj_zx("traj_zx", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig_traj_zx.set_window_properties({1150, 50}, {500, 500});
    fig_traj_zx.set_axes(IntervalVector({traj_z_bounds,traj_x_bounds}));

    Figure2D fig_xy("cut_xy", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig_xy.set_window_properties({50, 600}, {500, 500});
    fig_xy.set_axes(IntervalVector({x_bounds,y_bounds}));

    Figure2D fig_yz("cut_yz", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig_yz.set_window_properties({600, 600}, {500, 500});
    fig_yz.set_axes(IntervalVector({y_bounds,z_bounds}));

    Figure2D fig_zx("cut_zx", GraphicOutput::VIBES | GraphicOutput::IPE);
    fig_zx.set_window_properties({1150, 600}, {500, 500});
    fig_zx.set_axes(IntervalVector({z_bounds,x_bounds}));


    for (const auto &p : state_traj)
        if (p.first <= tf_discrete)
        {
            cout << "t: " << p.first  << endl;
            output.draw_sphere(p.second.subvector(0,2), Matrix::Identity(3,3)*0.025,Color::dark_green());
            fig_traj_xy.draw_circle({p.second[0], p.second[1]}, 0.02, StyleProperties({Color::dark_green(),Color::dark_green(0.5)},"z:2"));
            fig_traj_yz.draw_circle({p.second[1], p.second[2]}, 0.01, StyleProperties({Color::dark_green(),Color::dark_green(0.5)},"z:2"));
            fig_traj_zx.draw_circle({p.second[2], p.second[0]}, 0.005, StyleProperties({Color::dark_green(),Color::dark_green(0.5)},"z:2"));
        }
    fig_xy.draw_point({state_traj(tf_discrete)[0], state_traj(tf_discrete)[1]},StyleProperties({Color::dark_green(),Color::dark_green(0.5)},"z:2"));
    fig_yz.draw_point({state_traj(tf_discrete)[1], state_traj(tf_discrete)[2]},StyleProperties({Color::dark_green(),Color::dark_green(0.5)},"z:2"));
    fig_zx.draw_point({state_traj(tf_discrete)[2], state_traj(tf_discrete)[0]},StyleProperties({Color::dark_green(),Color::dark_green(0.5)},"z:2"));

    set_nb_threads(max_threads());
    auto wrappers = readIMapWrappers("../wrappers/dubins_wrappers.txt");


    double rad = 0.15;
    double x1_init = -5.0, x2_init = -0.45, theta_init=0.05;

    double theta_cut = 0.8;

    // Base

    VectorVar X(2);
    AnalyticFunction psi0 ({X},{rad*1/sqrt(1+sqr(X[0])+sqr(X[1])),rad*X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),rad*X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

    OctaSym id ({1,2,3});
    OctaSym s1 ({-2,1,3});
    OctaSym s2 ({3,2,-1});

    capd::IMap vectorField_wrap(std::get<0>(wrappers));

    auto peibos_output_base = PEIBOS(vectorField_wrap, tf_discrete, dt_discrete, psi0, {id,s1,s1*s1,s1.invert(),s2,s2.invert()}, 0.1, {x1_init, x2_init,theta_init}, true);

    auto m_v_par_base = reach_set(peibos_output_base);

    for (auto& [t, v_par] : m_v_par_base)
    {
        for (const auto& par: v_par)
            {
                output.draw_parallelepiped(par, StyleProperties(Color::pink(0.5),"base"));
                if (t == tf_discrete)
                    {
                        if (par.box()[2].intersects(Interval(state_traj(tf_discrete)[2])))
                            fig_xy.draw_zonotope(par.proj({0,1}), StyleProperties({Color::pink(), Color::pink(0.5)},"base","z:1"));
                        if (par.box()[0].intersects(Interval(state_traj(tf_discrete)[0])))
                            fig_yz.draw_zonotope(par.proj({1,2}), StyleProperties({Color::pink(), Color::pink(0.5)},"base","z:1"));
                        if (par.box()[1].intersects(Interval(state_traj(tf_discrete)[1])))
                            fig_zx.draw_zonotope(par.proj({2,0}), StyleProperties({Color::pink(), Color::pink(0.5)},"base","z:1"));
                    }
                fig_traj_xy.draw_zonotope(par.proj({0,1}), StyleProperties(Color::pink(0.5),"base_traj_xy"));
                fig_traj_yz.draw_zonotope(par.proj({1,2}), StyleProperties(Color::pink(0.5),"base_traj_yz"));
                fig_traj_zx.draw_zonotope(par.proj({2,0}), StyleProperties(Color::pink(0.5),"base_traj_zx"));
            }
    }
        

    // Reachability

    OctaSym id_reach ({1,2,3,4,5,6});
    OctaSym s1_reach ({-2,1,3,-5,4,6});
    OctaSym s2_reach ({3,2,-1,6,5,-4});

    AnalyticFunction psi0_reach ({X},{rad*1/sqrt(1+sqr(X[0])+sqr(X[1])),rad*X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),rad*X[1]/sqrt(1+sqr(X[0])+sqr(X[1])),1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

    capd::IMap vectorField_wrap_reach(std::get<2>(wrappers));

    // auto peibos_output_reach = PEIBOS(vectorField_wrap_reach,5,dt_discrete, psi0_reach, {id_reach, s1_reach, s1_reach * s1_reach, s1_reach.invert(), s2_reach, s2_reach.invert()} , 0.01, {x1_init,x2_init,theta_init,0,0,0}, true);
    vector<OctaSym> sym_vec_reach ({id_reach, s1_reach, s1_reach * s1_reach, s1_reach.invert(), s2_reach, s2_reach.invert()});

    int idx_sym = 0;
    for (auto& sym: sym_vec_reach)
    {
        double eps = 0.0125;
        auto peibos_output_reach = PEIBOS(vectorField_wrap_reach,tf_discrete,2.0, psi0_reach, {sym} , eps, {x1_init,x2_init,theta_init,0,0,0}, true);

        auto m_v_par_reach = reach_set(peibos_output_reach);

        for (auto& [t, v_par] : m_v_par_reach)
        {
            output.draw_sphere(state_traj.subvector(0, 2)(t), Matrix::Identity(3,3)*0.01 , StyleProperties({Color::green()},"trajectory"));
            for (const auto& par: v_par)
            {
                // output.draw_zonotope(par.proj({0,1,2}), StyleProperties(Color::red(0.5),"reach_index_"+to_string(idx_sym)));
                if (t == tf_discrete)
                {
                    if (par.box()[2].intersects(Interval(state_traj(tf_discrete)[2])))
                        fig_xy.draw_zonotope(par.proj({0,1}), StyleProperties({Color::red(),Color::red()},"reach_at_cut"));
                    if (par.box()[0].intersects(Interval(state_traj(tf_discrete)[0])))
                        fig_yz.draw_zonotope(par.proj({1,2}), StyleProperties({Color::red(),Color::red()},"reach_at_cut"));
                    if (par.box()[1].intersects(Interval(state_traj(tf_discrete)[1])))
                        fig_zx.draw_zonotope(par.proj({2,0}), StyleProperties({Color::red(),Color::red()},"reach_at_cut"));
                }
                    
            }
            
        }
        idx_sym++;
    }
}

#include <codac>
#include <codac-unsupported.h>
using namespace codac2;
using namespace std;

#include "cgal_tools.h"

Figure2D fig_inversion_in_out("inversion_in_out", GraphicOutput::VIBES | GraphicOutput::IPE);


static ColorMap peibos_cmap()
{
    ColorMap cmap(Model::HSV);
    int i = 0;
    for (int h = 300; h > 0; h -= 10)
    {
        cmap[i] = Color({(float)h, 50., 100., 50.}, Model::HSV);
        i++;
    }
    return cmap;
}

void triangulate(const vector<Parallelepiped> &v_par, const IntervalVector &init_box, const string &name)
{
    Figure2D fig_2d(name, GraphicOutput::VIBES | GraphicOutput::IPE);
    Figure2D fig_inside(name + "_polygon_inside", GraphicOutput::VIBES | GraphicOutput::IPE);
    Figure2D fig_outside(name + "_polygon_outside", GraphicOutput::VIBES | GraphicOutput::IPE);

    fig_2d.set_axes(axis(0, {init_box[0].lb(), init_box[0].ub()}), axis(1, {init_box[1].lb(), init_box[1].ub()}));
    fig_2d.set_window_properties({50, 50}, {500, 500});

    IntervalVector bbox(init_box);
    bbox.inflate(0.2);

    fig_inside.set_axes(bbox);
    fig_inside.set_window_properties({50, 600}, {500, 500});

    fig_outside.set_axes(bbox);
    fig_outside.set_window_properties({600, 600}, {500, 500});

    fig_outside.draw_box(init_box, StyleProperties(Color::red(), "w:0.02"));

    for (auto &p : v_par)
        fig_2d.draw_parallelepiped({p.c, p.A}, StyleProperties::boundary());

    Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_par);

    Polygon_2 outer_boundary = oriented_polygon(poly_with_holes.outer_boundary());
    fig_inside.draw_polygon(to_codac(outer_boundary), StyleProperties(Color::red(), "w:0.02"));
    fig_outside.draw_polygon(to_codac(outer_boundary), {Color::black(), Color::black(0.5)});

    for (auto it = poly_with_holes.holes_begin(); it != poly_with_holes.holes_end(); ++it)
    {
        Polygon_2 hole = oriented_polygon(*it);
        fig_inside.draw_polygon(to_codac(hole), {Color::black(), Color::black(0.5)});
    }

    vector<vector<Polygon>> inside_polygons;
    vector<vector<Polygon>> outside_polygons;

    vector<vector<Polygon>> holes_before = triangulate_holes(v_par);

    for (auto &hole : holes_before)
    {
        if (is_inside(hole, v_par))
        {
            inside_polygons.push_back(hole);
            for (auto &p : hole)
                fig_2d.draw_polygon(p, StyleProperties::inside());
        }

        else
        {
            outside_polygons.push_back(hole);
            for (auto &p : hole)
                fig_2d.draw_polygon(p, StyleProperties::outside());
        }
    }

    vector<Polygon> exterior_polygons = triangulate_outside(v_par, init_box);

    outside_polygons.push_back(exterior_polygons);

    for (auto &p : exterior_polygons)
        fig_2d.draw_polygon(p, StyleProperties::outside());
}

void triangulate(const vector<IntervalVector> &v_b, const IntervalVector &init_box)
{

    fig_inversion_in_out.draw_circle({0.5,0},1.75,Color::blue());

    fig_inversion_in_out.set_axes(axis(0, {init_box[0].lb(), init_box[0].ub()}), axis(1, {init_box[1].lb(), init_box[1].ub()}));
    fig_inversion_in_out.set_window_properties({1150, 50}, {500, 500});

    IntervalVector bbox(init_box);
    bbox.inflate(0.2);

    Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_b);

    Polygon_2 outer_boundary = oriented_polygon(poly_with_holes.outer_boundary());
    fig_inversion_in_out.draw_polygon(to_codac(outer_boundary), StyleProperties::boundary());
    cout << "Outer area: " << outer_boundary.area() << endl;
    for (auto it = poly_with_holes.holes_begin(); it != poly_with_holes.holes_end(); ++it)
    {
        Polygon_2 hole = oriented_polygon(*it);
        fig_inversion_in_out.draw_polygon(to_codac(hole), StyleProperties::inside());
        cout << "Hole area: " << hole.area() << endl;
    }
}

void triangulate_and_eliminate_fake(const vector<Parallelepiped> &v_par, const IntervalVector &init_box, const string &name)
{
    Figure2D fig_2d(name + "_fake", GraphicOutput::VIBES | GraphicOutput::IPE);
    Figure2D fig_no_fake(name + "_no_fake", GraphicOutput::VIBES | GraphicOutput::IPE);

    fig_2d.set_axes(axis(0, {init_box[0].lb(), init_box[0].ub()}), axis(1, {init_box[1].lb(), init_box[1].ub()}));
    fig_2d.set_window_properties({100, 100}, {500, 500});

    fig_no_fake.set_axes(axis(0, {init_box[0].lb(), init_box[0].ub()}), axis(1, {init_box[1].lb(), init_box[1].ub()}));
    fig_no_fake.set_window_properties({750, 100}, {500, 500});

    for (auto &p : v_par)
        fig_2d.draw_parallelepiped({p.c, p.A}, StyleProperties::boundary());

    vector<vector<Polygon>> inside_polygons;
    vector<vector<Polygon>> outside_polygons;

    vector<vector<Polygon>> holes_before = triangulate_holes(v_par);

    for (auto &hole : holes_before)
    {
        if (is_inside(hole, v_par))
        {
            inside_polygons.push_back(hole);
            for (auto &p : hole)
                fig_2d.draw_polygon(p, StyleProperties::inside());
        }

        else
        {
            outside_polygons.push_back(hole);
            for (auto &p : hole)
                fig_2d.draw_polygon(p, StyleProperties::outside());
        }
    }

    vector<Polygon> exterior_polygons = triangulate_outside(v_par, init_box);

    outside_polygons.push_back(exterior_polygons);

    for (auto &p : exterior_polygons)
    {
        fig_2d.draw_polygon(p, StyleProperties::outside());
    }

    auto v_par_no_fake = eliminate_fake_parallelepipeds(v_par, outside_polygons, inside_polygons);

    for (auto &p : v_par_no_fake)
        fig_no_fake.draw_parallelepiped({p.c, p.A}, StyleProperties::boundary());

    for (auto &v_poly : inside_polygons)
        for (auto &p : v_poly)
            fig_no_fake.draw_polygon(p, StyleProperties::inside());

    for (auto &v_poly : outside_polygons)
        for (auto &p : v_poly)
            fig_no_fake.draw_polygon(p, StyleProperties::outside());
}

void triangulate_with_diff(const vector<Parallelepiped> &v_par1, const vector<Parallelepiped> &v_par2, const IntervalVector &init_box, const string &name)
{
    Figure2D fig_2d(name, GraphicOutput::VIBES | GraphicOutput::IPE);
    Figure2D fig_poly_1(name + "_polygon_1", GraphicOutput::VIBES | GraphicOutput::IPE);
    Figure2D fig_poly_2(name + "_polygon_2", GraphicOutput::VIBES | GraphicOutput::IPE);

    fig_2d.draw_circle({0.5,0},1.75,Color::blue());

    fig_2d.set_axes(axis(0, {init_box[0].lb(), init_box[0].ub()}), axis(1, {init_box[1].lb(), init_box[1].ub()}));
    fig_2d.set_window_properties({1150, 600}, {500, 500});

    IntervalVector bbox(init_box);
    bbox.inflate(0.2);

    fig_poly_1.set_axes(bbox);
    fig_poly_1.set_window_properties({600, 600}, {500, 500});

    fig_poly_2.set_axes(bbox);
    fig_poly_2.set_window_properties({50, 600}, {500, 500});


    Polygon_with_holes_2 poly_with_holes_1 = generate_polygon_with_hole(v_par1);
    Polygon_with_holes_2 poly_with_holes_2 = generate_polygon_with_hole(v_par2);

    fig_poly_1.draw_polygon(to_codac(poly_with_holes_1.outer_boundary()), StyleProperties::boundary());


    auto hole_in = *poly_with_holes_2.holes_begin();
    auto polygon_out = poly_with_holes_2.outer_boundary();
    fig_poly_2.draw_polygon(to_codac(hole_in), StyleProperties::inside());
    fig_poly_2.draw_polygon(to_codac(polygon_out), StyleProperties::boundary());

    std::list<Polygon_with_holes_2> result_out;
    CGAL::difference(poly_with_holes_1.outer_boundary(), hole_in, std::back_inserter(result_out));
    for (const auto &poly_with_holes : result_out)
    {
        Polygon_2 outer_boundary = oriented_polygon(poly_with_holes.outer_boundary());
        fig_2d.draw_polygon(to_codac(outer_boundary), StyleProperties::boundary());
        cout << "Outer area: " << outer_boundary.area() << endl;
    }

    for (auto it = poly_with_holes_1.holes_begin(); it != poly_with_holes_1.holes_end(); ++it)
    {
        Polygon_2 hole = oriented_polygon(*it);
        std::list<Polygon_with_holes_2> result;
        CGAL::difference(hole, polygon_out, std::back_inserter(result));
        fig_poly_1.draw_polygon(to_codac(hole), StyleProperties::inside());
        for (const auto &poly_with_holes : result)
        {
            Polygon_2 outer_boundary = oriented_polygon(poly_with_holes.outer_boundary());
            fig_2d.draw_polygon(to_codac(outer_boundary), StyleProperties::inside());
            cout << "Hole area: " << outer_boundary.area() << endl;
        }
    }
}
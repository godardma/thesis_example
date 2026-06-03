#include <codac>
#include <codac-unsupported.h>
using namespace codac2;
using namespace std;

#include "cgal_tools.h"

static ColorMap peibos_cmap()
{
  ColorMap cmap( Model::HSV );
  int i = 0;
        for(int h = 300 ; h > 0 ; h-=10)
        {
          cmap[i]=Color({(float)h,50.,100.,50.},Model::HSV);
          i++;
        }
  return cmap;
}

void triangulate(const vector<Parallelepiped>& v_par, const IntervalVector& init_box, const string& name)
{
  Figure2D fig_2d(name, GraphicOutput::VIBES | GraphicOutput::IPE); 
  Figure2D fig_inside(name+"_polygon_inside", GraphicOutput::VIBES | GraphicOutput::IPE); 
  Figure2D fig_outside(name+"_polygon_outside", GraphicOutput::VIBES | GraphicOutput::IPE); 


  fig_2d.set_axes(axis(0,{init_box[0].lb(),init_box[0].ub()}),axis(1,{init_box[1].lb(),init_box[1].ub()}));
  fig_2d.set_window_properties({50,50},{500,500});

  IntervalVector bbox (init_box);
  bbox.inflate(0.2);

  fig_inside.set_axes(bbox);
  fig_inside.set_window_properties({50,600},{500,500});

  fig_outside.set_axes(bbox);
  fig_outside.set_window_properties({600,600},{500,500});

  fig_outside.draw_box(init_box, StyleProperties(Color::red(),"w:0.02"));

  for (auto& p : v_par)
      fig_2d.draw_parallelepiped({p.c, p.A}, StyleProperties::boundary());

    Polygon_with_holes_2 poly_with_holes = generate_polygon_with_hole(v_par);

    Polygon_2 outer_boundary = oriented_polygon(poly_with_holes.outer_boundary());
    fig_inside.draw_polygon(to_codac(outer_boundary), StyleProperties(Color::red(),"w:0.02"));
    fig_outside.draw_polygon(to_codac(outer_boundary), {Color::black(),Color::black(0.5)});

    for (auto it = poly_with_holes.holes_begin(); it != poly_with_holes.holes_end(); ++it) 
    {
        Polygon_2 hole = oriented_polygon(*it);
        fig_inside.draw_polygon(to_codac(hole), {Color::black(),Color::black(0.5)});
    }

  vector<vector<Polygon>> inside_polygons;
  vector<vector<Polygon>> outside_polygons;

  vector<vector<Polygon>> holes_before = triangulate_holes(v_par);

  for (auto& hole : holes_before)
  {
      if (is_inside(hole, v_par))
      {
          inside_polygons.push_back(hole);
          for (auto& p : hole)
              fig_2d.draw_polygon(p,StyleProperties::inside());
      }
          
      else
      {
          outside_polygons.push_back(hole);
          for (auto& p : hole)
              fig_2d.draw_polygon(p, StyleProperties::outside());
      }
          
  }
  
  vector<Polygon> exterior_polygons = triangulate_outside(v_par,init_box);

  outside_polygons.push_back(exterior_polygons);

  for (auto& p : exterior_polygons)
      fig_2d.draw_polygon(p, StyleProperties::outside());
}

void triangulate_and_eliminate_fake(const vector<Parallelepiped>& v_par, const IntervalVector& init_box, const string& name)
{
  Figure2D fig_2d(name + "_fake", GraphicOutput::VIBES | GraphicOutput::IPE); 
  Figure2D fig_no_fake(name + "_no_fake", GraphicOutput::VIBES | GraphicOutput::IPE);

  fig_2d.set_axes(axis(0,{init_box[0].lb(),init_box[0].ub()}),axis(1,{init_box[1].lb(),init_box[1].ub()}));
  fig_2d.set_window_properties({100,100},{500,500});

  fig_no_fake.set_axes(axis(0,{init_box[0].lb(),init_box[0].ub()}),axis(1,{init_box[1].lb(),init_box[1].ub()}));
  fig_no_fake.set_window_properties({750,100},{500,500});

  for (auto& p : v_par)
      fig_2d.draw_parallelepiped({p.c, p.A}, StyleProperties::boundary());

  vector<vector<Polygon>> inside_polygons;
  vector<vector<Polygon>> outside_polygons;

  vector<vector<Polygon>> holes_before = triangulate_holes(v_par);

  for (auto& hole : holes_before)
  {
      if (is_inside(hole, v_par))
      {
          inside_polygons.push_back(hole);
          for (auto& p : hole)
              fig_2d.draw_polygon(p,StyleProperties::inside());
      }
          
      else
      {
          outside_polygons.push_back(hole);
          for (auto& p : hole)
              fig_2d.draw_polygon(p, StyleProperties::outside());
      }
          
  }
  
  vector<Polygon> exterior_polygons = triangulate_outside(v_par,init_box);

  outside_polygons.push_back(exterior_polygons);

  for (auto& p : exterior_polygons)
  {
      fig_2d.draw_polygon(p, StyleProperties::outside());
  }

  auto v_par_no_fake = eliminate_fake_parallelepipeds(v_par, outside_polygons, inside_polygons);

  for (auto& p : v_par_no_fake)
      fig_no_fake.draw_parallelepiped({p.c, p.A}, StyleProperties::boundary());

  for (auto& v_poly : inside_polygons)
      for (auto& p : v_poly)
          fig_no_fake.draw_polygon(p, StyleProperties::inside());

  for (auto& v_poly : outside_polygons)
      for (auto& p : v_poly)
          fig_no_fake.draw_polygon(p, StyleProperties::outside());
}
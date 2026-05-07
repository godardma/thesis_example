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

  fig_2d.set_axes(axis(0,{init_box[0].lb(),init_box[0].ub()}),axis(1,{init_box[1].lb(),init_box[1].ub()}));
  fig_2d.set_window_properties({50,50},{500,500});

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
      fig_2d.draw_polygon(p, StyleProperties::outside());
}

void triangulate_and_eliminate_fake(const vector<Parallelepiped>& v_par, const IntervalVector& init_box, const string& name)
{
  Figure2D fig_2d(name + " Fake", GraphicOutput::VIBES | GraphicOutput::IPE); 
  Figure2D fig_no_fake(name + " No Fake", GraphicOutput::VIBES | GraphicOutput::IPE);

  fig_2d.set_axes(axis(0,{init_box[0].lb(),init_box[0].ub()}),axis(1,{init_box[1].lb(),init_box[1].ub()}));
  fig_2d.set_window_properties({100,100},{500,500});

  fig_no_fake.set_axes(axis(0,{init_box[0].lb(),init_box[0].ub()}),axis(1,{init_box[1].lb(),init_box[1].ub()}));
  fig_no_fake.set_window_properties({750,100},{500,500});

  for (auto& p : v_par)
      fig_2d.draw_parallelepiped({p.c, p.A}, StyleProperties({Color::yellow(), Color::yellow(0.5)},"polygons"));

  vector<vector<Polygon>> inside_polygons;
  vector<vector<Polygon>> outside_polygons;

  vector<vector<Polygon>> holes_before = triangulate_holes(v_par);

  for (auto& hole : holes_before)
  {
      if (is_inside(hole, v_par))
      {
          inside_polygons.push_back(hole);
          for (auto& p : hole)
              fig_2d.draw_polygon(p,StyleProperties({Color::green(), Color::green(0.5)},"inside"));
      }
          
      else
      {
          outside_polygons.push_back(hole);
          for (auto& p : hole)
              fig_2d.draw_polygon(p, StyleProperties({Color::red(), Color::purple(0.5)},"outside"));
      }
          
  }
  
  vector<Polygon> exterior_polygons = triangulate_outside(v_par,init_box);

  outside_polygons.push_back(exterior_polygons);

  for (auto& p : exterior_polygons)
  {
      fig_2d.draw_polygon(p, {Color::red(), Color::red(0.5)});
  }

  auto v_par_no_fake = eliminate_fake_parallelepipeds(v_par, outside_polygons, inside_polygons);

  for (auto& p : v_par_no_fake)
      fig_no_fake.draw_parallelepiped({p.c, p.A}, StyleProperties({Color::yellow(), Color::yellow(0.5)},"polygons"));

  for (auto& v_poly : inside_polygons)
      for (auto& p : v_poly)
          fig_no_fake.draw_polygon(p, StyleProperties({Color::green(), Color::green(0.5)},"inside"));

  for (auto& v_poly : outside_polygons)
      for (auto& p : v_poly)
          fig_no_fake.draw_polygon(p, StyleProperties({Color::red(), Color::red(0.5)},"outside"));
}
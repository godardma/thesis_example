// Author: Maël GODARD

#include <codac>
#include <codac-unsupported.h>
#include <codac-capd.h>

#include <chrono>

using namespace std;
using namespace codac2;

ColorMap rainbow_05(float alpha=1.)
{
  ColorMap cmap( Model::HSV );
  int i = 0;
  for(int h = 300 ; h > 0 ; h-=10)
  {
    cmap[i]=Color({(float)h,50.,100.,100.f*alpha},Model::HSV);
    i++;
  }
  return cmap;
}

bool verify_criteria (Parallelepiped& p)
{
  IntervalVector bbox = p.box();
  CtcDist c;
  Interval d(0.03,+oo), a1(1), a2(0);
  c.contract(bbox[0], bbox[1], a1, a2, d);
  return bbox.is_empty();

}

  Parallelepiped PEIBOS3D_step(const capd::IMap& IMap, double tf, IntervalVector X, const AnalyticFunction<VectorType>& psi_0, const OctaSym& symmetry, double epsilon, const Vector& offset)
  {
    int m = psi_0.input_size();
    int n = psi_0.output_size();

    assert(offset.size() == n);
    assert (m < n);
    assert (symmetries.size() > 0 && (int) symmetries[0].size() == n);

    clock_t t_start = clock();

    vector<Parallelepiped> output;
    
    // CAPD solver setup
    capd::IMap g (IMap);
    capd::IOdeSolver solver(g, 30);
    
    solver.setAbsoluteTolerance(1e-20);
    solver.setRelativeTolerance(1e-20);

    capd::ITimeMap timeMap(solver);
    capd::ITimeMap timeMap_punc(solver);

    capd::interval initialTime(0.);
    capd::interval finalTime(tf);


    // To get the flow function and its Jacobian (monodromy matrix) for [x]
    IntervalVector Y = symmetry(psi_0.eval(X)) + offset;

    capd::IMatrix monodromyMatrix(n,n);
    capd::ITimeMap::SolutionCurve solution(initialTime); 
    capd::IVector c =to_capd(Y);

    capd::C1Rect2Set s(c);
    timeMap(finalTime, s, solution);
    capd::IVector result = timeMap(finalTime, s, monodromyMatrix);
    auto JJf=to_codac(monodromyMatrix);

    // To get the flow function and its Jacobian (monodromy matrix) for x_hat
    auto xc = X.mid();
    auto yc = (symmetry(psi_0.eval(xc)) + offset).mid();

    capd::IMatrix monodromyMatrix_punc(n,n);
    capd::ITimeMap::SolutionCurve solution_punct(initialTime);
    capd::IVector c_punct =to_capd(IntervalVector(yc));

    capd::C1Rect2Set s_punct(c_punct);
    timeMap_punc(finalTime, s_punct, solution_punct);      
    capd::IVector result_punct = timeMap_punc(finalTime, s_punct, monodromyMatrix_punc);
    auto JJf_point=to_codac(monodromyMatrix_punc);

    // Center of the parallelepiped
    Vector z = Vector(to_codac(result).mid());
    
    auto p = parallelepiped_inclusion(to_codac(result_punct), JJf, JJf_point.mid(), psi_0, symmetry, X);

    return p;
  }

template <typename T>
void PEIBOS_recursive (const capd::IMap& i_map_wrapper, double tf, AnalyticFunction<T>& psi_0, OctaSym symmetry , IntervalVector X, double epsilon, Vector offset, Figure3D& figure_3d, vector<Parallelepiped>& output, double& epsilon_min)
{
  if (epsilon < 0.01)
  {
    cout << "The epsilon limit has been reached" << endl;
    exit(0);
  }

  epsilon_min = min(epsilon_min, epsilon);
  
  for (double t1 = X[0].lb(); t1 < X[0].ub(); t1 += epsilon)
  {
    for (double t2 = X[1].lb();t2 < X[1].ub(); t2+=epsilon)
    {
      IntervalVector X0({{t1,t1+epsilon},{t2,t2+epsilon}});
      auto parallel = PEIBOS3D_step(i_map_wrapper, tf, X0, psi_0, symmetry, epsilon, offset);
      if (verify_criteria(parallel))
      {
        output.push_back(parallel);
        Zonotope z = parallel.proj({0,1,2});

        Matrix A_proj (3,3);
        A_proj.col(0) = z.A.col(0);
        A_proj.col(1) = z.A.col(1);
        A_proj.col(2) = z.A.col(2);
        figure_3d.draw_parallelepiped(Parallelepiped(z.c, A_proj), rainbow_05(0.5).color(epsilon/0.5));
      }
      else
      {
        PEIBOS_recursive(i_map_wrapper, tf, psi_0, symmetry, X0, epsilon/2., offset, figure_3d, output, epsilon_min);
      }
    
    }
  }
}

template <typename T>
void PEIBOS_adaptative(const capd::IMap& i_map_wrapper, double tf, AnalyticFunction<T>& psi_0,const vector<OctaSym>& symmetries , double epsilon, Vector offset, Figure3D& figure_3d, vector<Parallelepiped>& output)
{
  double epsilon_init = epsilon;
  ColorMap cmap = rainbow_05(0.5);

    for (int i = 0; i < symmetries.size(); i++)
    {
      IntervalVector X = IntervalVector({{-1,1},{-1,1}});
      PEIBOS_recursive(i_map_wrapper, tf, psi_0, symmetries[i], X, epsilon, offset, figure_3d,output, epsilon_init);
    }
  cout << "Minimum epsilon reached: " << epsilon_init << endl;
}

int main()
{


  
  capd::IMap vectorField_wrap("par:v,T;var:x1,x2,x3,t;fun:v*cos(x3),v*sin(x3),cos(2*3.14*t/T),1;");
  vectorField_wrap.setParameter("v",1.0);
  vectorField_wrap.setParameter("T",1.0);
  
  VectorVar X_3d(2);
  AnalyticFunction psi0 ({X_3d},{0.02*1/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),0.02*X_3d[0]/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),0.02*X_3d[1]/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),0});
  OctaSym id_3d ({1,2,3,4});
  OctaSym s1 ({-2,1,3,4});
  OctaSym s2 ({3,2,-1,4});

  double epsilon = 0.5;

  Figure3D output ("dubins_reach");
  output.draw_axes(0.05,{1,-0.05,-0.05});

  Figure3D output_traj ("dubins_reach_traj");
  output_traj.draw_axes(0.2,{-0.05,-0.05,-0.05});

  Figure2D fig ("project_dubins",GraphicOutput::VIBES|GraphicOutput::IPE);
  fig.set_window_properties({50,50},{500,500});
  fig.set_axes(axis(0,{0.95,1.05}), axis(1,{-0.05,0.05}));

  fig.draw_ring({1,0},{0.03,0.05}, {Color({255,102,0,255}),Color({255,102,0,130})});

  Figure2D fig_fine ("project_dubins_fine",GraphicOutput::VIBES|GraphicOutput::IPE);
  fig_fine.set_window_properties({650,50},{500,500});
  fig_fine.set_axes(axis(0,{0.95,1.05}), axis(1,{-0.05,0.05}));

  fig_fine.draw_ring({1,0},{0.03,0.05}, {Color({255,102,0,255}),Color({255,102,0,130})});

  vector<Parallelepiped> v_par ;

  // start the timer
  auto start = chrono::high_resolution_clock::now();
  PEIBOS_adaptative(vectorField_wrap, 1.0, psi0, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, epsilon, {0.0,0.0,0.0,0.0}, output, v_par);

  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
  cout << "Execution time  with adpative: " << duration.count()/1000. << " seconds" << endl;

  for (const auto& p : v_par)
    fig.draw_zonotope(p.proj({0,1}), {Color::black(),Color::green(0.5)});

  auto peibos_output_traj = PEIBOS(vectorField_wrap, 1.0, 0.2, psi0, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, 0.25, {0.0,0.0,0.0,0.0});
  auto m_v_par_traj = reach_set(peibos_output_traj);

  for (auto& [t, v_par] : m_v_par_traj)
    for (const auto& par: v_par)
      output_traj.draw_zonotope(par.proj({0,1,2}), {Color::blue(0.5)});

  start = chrono::high_resolution_clock::now();
  auto peibos_output_fine = PEIBOS(vectorField_wrap, 1.0, 1.0, psi0, {id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()}, 0.03125, {0.0,0.0,0.0,0.0});
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
  cout << "Execution time without adpative: " << duration.count()/1000. << " seconds" << endl;
  auto m_v_par_fine = reach_set(peibos_output_fine);

  for (auto& [t, v_par] : m_v_par_fine)
    if (t == 1.0)
      for (const auto& par: v_par)
        fig_fine.draw_zonotope(par.proj({0,1}), {Color::black(),Color::green(0.5)});

}

  

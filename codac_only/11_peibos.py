from codac import *
import numpy as np

def subdivide(x, epsilon, Lx):
  if x.max_diam()<=epsilon:
    Lx.append(x)
  else:
    (x1,x2)=x.bisect_largest(0.5)
    subdivide(x1,epsilon,Lx)
    subdivide(x2,epsilon,Lx)

def box_enclosure(f,psi_0,Sigma,epsilon):
  Lb = []
  X = VectorVar(2)
  for sigma in Sigma:
    g_i = AnalyticFunction([X],f(sigma(psi_0(X))))
    Lx = []
    subdivide(IntervalVector.constant(2,Interval(-1,1)),epsilon,Lx)
    for x in Lx:
      b = g_i.eval(EvalMode.NATURAL,x)
      Lb.append(b)
  return Lb

set_nb_threads(8)

# 3D example of the PEIBOS algorithm

y_3d = VectorVar(3)
f_3d = AnalyticFunction([y_3d],[sqr(y_3d[0])-sqr(y_3d[1])+y_3d[0],2*y_3d[0]*y_3d[1]+y_3d[1],y_3d[2]])

X_3d = VectorVar(2)
psi0_3d = AnalyticFunction([X_3d],[1/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),X_3d[0]/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),X_3d[1]/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1]))])

id_3d = OctaSym([1, 2, 3])
s1 = OctaSym([-2, 1, 3])
s2 = OctaSym([3, 2, -1])

# for computation time comparison
# epsilons = [0.0625,0.25,0.125,0.0625,0.03125,1/2**6,1/2**7,1/2**8]

# for display
epsilons = [0.0625,0.25,0.125,0.0625]

figure_3d = Figure3D("PEIBOS")
figure_3d.draw_axes(2.5)

for i in range (len(epsilons)):
  print(epsilons[i])
  v_par_3d = PEIBOS(f_3d,psi0_3d,[id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()],epsilons[i])
  v_box_3d = box_enclosure(f_3d,psi0_3d,[id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()],epsilons[i])
  index = 0
  volume_par = 0
  for p in v_par_3d:
    LUdec = IntvFullPivLU(p.A)
    volume_par += abs(LUdec.determinant())
    if i==0:
      # draw each sixth of the list in a different color for the first figure
      if index < len(v_par_3d)/6:
        color = Color.red(0.8)
      elif index < 2*len(v_par_3d)/6:
        color = Color.blue(0.8)
      elif index < 3*len(v_par_3d)/6:
        color = Color.red(0.8)
      elif index < 4*len(v_par_3d)/6:
        color = Color.blue(0.8)
      elif index < 5*len(v_par_3d)/6:
        color = Color.yellow(0.8)
      else:
        color = Color.yellow(0.8)
      figure_3d.draw_parallelepiped(p,StyleProperties(color,"atlas"))
      index += 1
    elif i<4:
      figure_3d.draw_parallelepiped(p,StyleProperties(Color.green(0.5),"epsilon_"+str(epsilons[i]).replace(".","_")))

  volume_box = 0
  for b in v_box_3d:
    volume_box += b.volume()
    if i==2:
      figure_3d.draw_box(b,StyleProperties(Color.blue(0.5),"boxes_epsilon_"+str(epsilons[i]).replace(".","_")))

  print(volume_par.mid())
  print(volume_box)
  print()

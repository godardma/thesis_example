from codac import *
import numpy as np

set_nb_threads(8)

# 3D example of the PEIBOS algorithm

y_3d = VectorVar(3)
f_3d = AnalyticFunction([y_3d],[sqr(y_3d[0])-sqr(y_3d[1])+y_3d[0],2*y_3d[0]*y_3d[1]+y_3d[1],y_3d[2]])

X_3d = VectorVar(2)
psi0_3d = AnalyticFunction([X_3d],[1/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),X_3d[0]/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1])),X_3d[1]/sqrt(1+sqr(X_3d[0])+sqr(X_3d[1]))])

id_3d = OctaSym([1, 2, 3])
s1 = OctaSym([-2, 1, 3])
s2 = OctaSym([3, 2, -1])


epsilons = [0.0625,0.25,0.125,0.0625,0.03125,1/2**6,1/2**7,1/2**8,1/2**9,1/2**10]

figure_3d = Figure3D("PEIBOS")
figure_3d.draw_axes(1.2)

for i in range (len(epsilons)):
  print(epsilons[i])
  v_par_3d = PEIBOS(f_3d,psi0_3d,[id_3d,s1,s1*s1,s1.invert(),s2,s2.invert()],epsilons[i])
  index = 0
  volume = 0
  for p in v_par_3d:
    LUdec = IntvFullPivLU(p.A)
    volume += abs(LUdec.determinant())
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

  print(volume.mid())
  print()

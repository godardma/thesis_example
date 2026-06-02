from codac import *

def split(x, epsilon, Lx):
  if x.max_diam()<=epsilon:
    Lx.append(x)
  else:
    (x1,x2)=x.bisect_largest(0.5)
    split(x1,epsilon,Lx)
    split(x2,epsilon,Lx)

def color(i):
  if i == 0:
    return Color.red()
  else :
    return Color.green(0.5)
  
def PEIB(f,x0,epsilon):
  Lp = []
  Lx = []
  split(x0,epsilon,Lx)
  for x in Lx:
    p = f.parallelepiped_eval(x)
    Lp.append(p)
    
  return Lx, Lp
  
  

X = VectorVar(2)
psi_0 = AnalyticFunction([X],[1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))])

x0 = IntervalVector.constant(2,Interval(-1,1))
epsilons = [0.02,0.25,0.125,0.0625]

fig0 = Figure2D("parallelepiped_enclosure_epsilon_"+str(epsilons[0]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig0.set_window_properties([50,50],[500,500])
fig0.set_axes(x0)

fig1 = Figure2D("parallelepiped_enclosure_epsilon_"+str(epsilons[1]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig1.set_window_properties([600,50],[500,500])
fig1.set_axes(x0)

fig2 = Figure2D("parallelepiped_enclosure_epsilon_"+str(epsilons[2]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig2.set_window_properties([50,600],[500,500])
fig2.set_axes(x0)

fig3 = Figure2D("parallelepiped_enclosure_epsilon_"+str(epsilons[3]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig3.set_window_properties([600,600],[500,500])
fig3.set_axes(x0)

figs = [fig0,fig1,fig2,fig3]
fig3d = Figure3D("parallelepiped_enclosure")
fig3d.draw_axes(0.8)

for i in range (len(epsilons)):
  Lx,Lp =PEIB(psi_0,x0,epsilons[i])
  print(i)
  for x in Lx:
    figs[i].draw_box(x,[Color.black(),color(i)])
    print(x.max_diam())
  for p in Lp:
    fig3d.draw_parallelepiped(p,StyleProperties(color(i),"epsilon_"+str(epsilons[i]).replace(".","_")))

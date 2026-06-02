from codac import *

def split(x, epsilon, Lx):
  if x.max_diam()<=epsilon:
    Lx.append(x)
  else:
    (x1,x2)=x.bisect_largest(0.5)
    split(x1,epsilon,Lx)
    split(x2,epsilon,Lx)
  
def PEIB(f,x0,epsilon):
  Lp = []
  Lx = []
  split(x0,epsilon,Lx)
  for x in Lx:
    p = f.parallelepiped_eval(x)
    Lp.append(p)
  return Lx, Lp
  
  

X = VectorVar(1)
f = AnalyticFunction([X],[X[0],Interval(1.1,1.2)*sqr(X[0])])

x0 = IntervalVector.constant(1,Interval(-2,2))
epsilons = [1.,0.5,0.25,0.125]

y0 = IntervalVector([[-3.5,3.5],[-1.5,6.5]])

fig0 = Figure2D("parallelepiped_enclosure_set_valued_epsilon_"+str(epsilons[0]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig0.set_window_properties([50,50],[500,500])
fig0.set_axes(y0)

fig1 = Figure2D("parallelepiped_enclosure_set_valued_epsilon_"+str(epsilons[1]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig1.set_window_properties([600,50],[500,500])
fig1.set_axes(y0)

fig2 = Figure2D("parallelepiped_enclosure_set_valued_epsilon_"+str(epsilons[2]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig2.set_window_properties([50,600],[500,500])
fig2.set_axes(y0)

fig3 = Figure2D("parallelepiped_enclosure_set_valued_epsilon_"+str(epsilons[3]).replace(".","_"), GraphicOutput.VIBES | GraphicOutput.IPE)
fig3.set_window_properties([600,600],[500,500])
fig3.set_axes(y0)

figs = [fig0,fig1,fig2,fig3]



for i in range (len(epsilons)):
  Lx,Lp =PEIB(f,x0,epsilons[i])
  for p in Lp:
    figs[i].draw_parallelepiped(p,[Color.black(),Color.green(0.5)])


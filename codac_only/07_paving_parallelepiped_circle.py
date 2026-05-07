from codac import *

x = VectorVar(1)
f = AnalyticFunction([x], [cos(x[0]),sin(x[0])])

L_p = []
nb_parallelepiped = 15
delta_angle = 2*PI/nb_parallelepiped

ctc_union = CtcUnion(2)

for i in range (nb_parallelepiped):
  X = IntervalVector([Interval(i*delta_angle,(i+1)*delta_angle)])
  p = f.parallelepiped_eval(X)
  L_p.append(p)
  ctc_union |= CtcParallelepiped(p)

x0 = IntervalVector([[-1.5,1.5],[-1.5,1.5]]) 
p = pave(x0,ctc_union,0.15)

# Graphical setup

fig1 = Figure2D("Paving_of_circle_parallelelepiped", GraphicOutput.VIBES | GraphicOutput.IPE)
fig1.set_window_properties([50,50],[500,500])
fig1.set_axes(x0)

fig2 = Figure2D("Paving_of_circle_parallelelepiped_in_out", GraphicOutput.VIBES | GraphicOutput.IPE)
fig2.set_window_properties([600,50],[500,500])
fig2.set_axes(x0)

fig1.draw_paving(p)

for par in L_p:
  fig1.draw_parallelepiped(par,StyleProperties(Color.red(),"w:0.008","--"))
  
v_cs = p.connected_subsets(PavingOut.outer_complem)

for i in range (len(v_cs)):
  if i == 0:
    fig2.draw_subpaving(v_cs[i],StyleProperties.outside())
  else:
    fig2.draw_subpaving(v_cs[i],StyleProperties.inside())

for box in p.boxes(PavingOut.outer):
  fig2.draw_box(box,StyleProperties.boundary())
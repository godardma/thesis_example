from codac import *

x = VectorVar(2)
f = AnalyticFunction([x],sqrt(sqr(x[0])+sqr(x[1])))

x0 = IntervalVector([[-4,4],[-4,4]])

epsilon = 0.1
p = sivia(x0,f,Interval(1,3),epsilon)

fig = Figure2D("SIVIA", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500])
fig.set_axes(x0)

fig.draw_paving(p)

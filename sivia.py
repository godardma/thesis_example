from codac import *

x = VectorVar(2)
f = AnalyticFunction([x],sqrt(sqr(x[0])+sqr(x[1])))

x0 = IntervalVector([[-4,4],[-4,4]])
p = sivia(x0,f,Interval(1,3),0.1)



fig = Figure2D("SIVIA", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500]) # set the window position and size
fig.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

fig.draw_paving(p)
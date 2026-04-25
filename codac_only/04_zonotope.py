from codac import *

z = Vector([4,1.5])
Az = Matrix([[-0.2,0.2  ,0.08],
             [0.1,0.04,0.18]])

fig = Figure2D("Zonotope", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500]) # set the window position and size
x0 = IntervalVector([[3,5],[1,2]])
fig.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1


fig.draw_zonotope(Zonotope(z,Az),StyleProperties([Color.black(),Color.gray(0.5)],"zonotope","w:0.02"))

for j in range (Az.cols()):
  if j == 0:
    fig.draw_arrow(z,z+Az.col(j),0.1,[Color.red(),Color.red()])
  if j == 1:
    fig.draw_arrow(z,z+Az.col(j),0.1,[Color.green(),Color.green()])
  if j == 2:
    fig.draw_arrow(z,z+Az.col(j),0.1,[Color.blue(),Color.blue()])
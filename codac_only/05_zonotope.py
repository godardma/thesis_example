from codac import *

z1 = Vector([4,1.5])
A1 = Matrix([[-0.2, 0.2, 0.08],
             [0.1, 0.04, 0.18]])
Z1 = Zonotope(z1,A1)

z2 = Vector([3,-0.5])
A2 = Matrix([[0.2, -0.1],
             [0.15,0.01]])
Z2 = Zonotope(z2,A2)

fig1 = Figure2D("Zonotope", GraphicOutput.VIBES | GraphicOutput.IPE)
fig1.set_window_properties([50,50],[500,500]) 
fig1.set_axes(IntervalVector([[3,5],[1,2]]))


fig1.draw_zonotope(Z1,StyleProperties([Color.black(),Color.gray(0.5)],"zonotope","w:0.02"))

for j in range (A1.cols()):
  if j == 0:
    fig1.draw_arrow(z1,z1+A1.col(j),0.1,[Color.red(),Color.red()])
  if j == 1:
    fig1.draw_arrow(z1,z1+A1.col(j),0.1,[Color.green(),Color.green()])
  if j == 2:
    fig1.draw_arrow(z1,z1+A1.col(j),0.1,[Color.blue(),Color.blue()])

fig2 = Figure2D("Zonotope Sum", GraphicOutput.VIBES | GraphicOutput.IPE)
fig2.set_window_properties([600,50],[500,500]) 
fig2.set_axes(IntervalVector([[1,9],[-1,2.3]]))


fig2.draw_zonotope(Z1,StyleProperties([Color.black(),Color.gray(0.5)],"zonotope","w:0.02"))
fig2.draw_zonotope(Z2,StyleProperties([Color.black(),Color.blue(0.5)],"zonotope","w:0.02"))
fig2.draw_zonotope(Z1+Z2,StyleProperties([Color.black(),Color.red(0.5)],"zonotope","w:0.02"))
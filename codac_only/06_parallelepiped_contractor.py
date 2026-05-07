from codac import *

def CtcParallelepiped(x):
  global p
  x_p = inverse_enclosure(p.A)*(x-p.c)
  x_p &= IntervalVector.constant(p.c.size(),Interval(-1,1))
  x_c = p.c + p.A*x_p
  return x_c & x


def pave_contract(ctc, x, epsilon, X_plus):
  x = ctc(x)
  if not (x.is_empty()):
    if x.max_diam()<epsilon:
      X_plus.append(x)
    else:
      (x1,x2)=x.bisect_largest()
      pave_contract(ctc, x1, epsilon, X_plus)
      pave_contract(ctc, x2, epsilon, X_plus)


p = Parallelepiped([1.5,2.8],Matrix([[0.5,0.4],[-0.1,0.2]]))
x0 = IntervalVector([[0,3],[2,3.5]])


X_plus_contract = []
pave_contract(CtcParallelepiped, IntervalVector(x0), 0.1, X_plus_contract)

# Graphical setup

fig1 = Figure2D("Paving_of_CtcParallelepiped", GraphicOutput.VIBES | GraphicOutput.IPE)
fig1.set_window_properties([50,50],[500,500])
fig1.set_axes(x0)

fig1.draw_box(x0,StyleProperties.outside())
fig1.draw_parallelepiped(p,StyleProperties([Color.red(),Color.red(0.2)],"w:0.01","--"))

area_contract = 0
for box in X_plus_contract:
  area_contract += box.volume()
  fig1.draw_box(box,StyleProperties.boundary())
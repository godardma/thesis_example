from codac import *

a = 11; b = 7
x = VectorVar(2)
f = AnalyticFunction([x], sqr(sqr(x[0])+x[1]-a)+sqr(x[0]+sqr(x[1])-b))

ctc_inv = CtcInverse (f, 80)

x0 = IntervalVector.constant(2,[-7,7])
p = pave(x0, ctc_inv,0.01)

boxes_to_contract = [IntervalVector([[-3,6.5],[-5,-2]]),
                     IntervalVector([[-3,-1],[0,1.5]]),
                     IntervalVector([[4,6],[0.5,2]])]

fig = Figure2D("Himmelblau", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500]) # set the window position and size
fig.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

# fig.draw_paving(p)

L_bound = p.boxes(PavingOut.outer)

for box in L_bound:
  fig.draw_box(box,[Color.black(),Color.black()])

for box in boxes_to_contract:
  fig.draw_box(box,Color.blue())
  box = ctc_inv.contract(box)
  fig.draw_box(box,Color.red())
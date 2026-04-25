from codac import *


a1,a2,b1,b2 = 2,3,4,5
x = VectorVar(2)
f = AnalyticFunction([x], sqrt(sqr(x[0]-a1)+sqr(x[1])-a2))

ctc_inv = CtcInverse (f, Interval(2,3))

x0 = IntervalVector([[-6,7],[-6,6]])
p = pave(x0, ctc_inv,0.1)

# boxes_to_contract = [IntervalVector([[-3,6.5],[-5,-2]]),
#                      IntervalVector([[-4,-2],[0,1.5]]),
#                      IntervalVector([[4,6],[0.5,2]]),
#                      IntervalVector([[-0.5,0.5],[-0.5,0.5]]),
#                      IntervalVector([[2.5,3.5],[0.5,2.5]])]

fig = Figure2D("constraint_programming", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500]) # set the window position and size
fig.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

fig.draw_paving(p)
# L_bound = p.boxes(PavingOut.outer)

# for box in L_bound:
#   fig.draw_box(box,[Color.black(),Color.black()])

# for box in boxes_to_contract:
#   fig.draw_box(box,StyleProperties(Color.blue(),"w:0.05"))
#   box = ctc_inv.contract(box)
#   fig.draw_box(box,StyleProperties(Color.red(),"w:0.05"))
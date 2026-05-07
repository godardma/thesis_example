from codac import *


a_1, a_2 = Vector([2,3]),Vector([-2,-2])
r1,theta1 = Interval(2,3),Interval(-PI,-PI/2)
r2,theta2 = Interval(3,4.1),Interval(0,PI/2)
x = VectorVar(2)
y1,y2 = x-a_1,x-a_2
f = AnalyticFunction([x], [sqrt(sqr(y1[0])+sqr(y1[1])),atan2(y1[1],y1[0]),sqrt(sqr(y2[0])+sqr(y2[1])),atan2(y2[1],y2[0])])

ctc_inv = CtcInverse (f, IntervalVector([r1,theta1,r2,theta2]),False)
box_to_contract = IntervalVector([[0.5,1.5],[1,2.5]])



fig = Figure2D("constraint_programming", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500]) 
x0 = IntervalVector([[-3,4],[-3,4]])
fig.set_axes(x0)

fig.draw_point(a_1)
fig.draw_point(a_2)
fig.draw_pie(a_1,r1|0,theta1, Color.gray())
fig.draw_pie(a_2,r2|0,theta2, Color.gray())
fig.draw_pie(a_1,r1,theta1,[Color.black(),Color.gray(0.15)])
fig.draw_pie(a_2,r2,theta2,[Color.black(),Color.gray(0.15)])


box_inv = IntervalVector(box_to_contract)
ctc_inv.contract(box_inv)

ctc_polar = CtcPolar()
box_pol1 = IntervalVector(box_to_contract)
box_pol2 = IntervalVector(box_to_contract)
box_pol1-=a_1
box_pol2-=a_2
ctc_polar.contract(box_pol1[0],box_pol1[1],r1,theta1)
ctc_polar.contract(box_pol2[0],box_pol2[1],r2,theta2)
box_pol1+=a_1
box_pol2+=a_2


fig.draw_box(box_to_contract,StyleProperties(Color.blue(),"w:0.05"))
fig.draw_box(box_inv,StyleProperties(Color.red(),"w:0.05"))
fig.draw_box(box_pol1&box_pol2,StyleProperties(Color.green(),"w:0.05",".."))

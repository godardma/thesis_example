from codac import *
import numpy as np

import time

def CtcTryskell(x):
  x_c = IntervalVector.empty(2)
  ctcPolar = CtcPolar()
  
  x_O_R2_1 = IntervalVector(x) + O
  ctcPolar.contract(x_O_R2_1[0], x_O_R2_1[1], Interval(R2), Interval(theta_O_R2_1))
  x_c |= x_O_R2_1 - O

  x_O_R2_2 = IntervalVector(x) + O
  ctcPolar.contract(x_O_R2_2[0], x_O_R2_2[1], Interval(R2), Interval(theta_O_R2_2))
  x_c |= x_O_R2_2 - O

  x_O_R2_3 = IntervalVector(x) + O
  ctcPolar.contract(x_O_R2_3[0], x_O_R2_3[1], Interval(R2), Interval(theta_O_R2_3))
  x_c |= x_O_R2_3 - O


  x_B1_r2_1 = IntervalVector(x) - B1
  ctcPolar.contract(x_B1_r2_1[0],x_B1_r2_1[1], Interval(r2), Interval(theta_B1_r2_1))
  x_c |= (x_B1_r2_1 + B1)

  x_B2_r2_1 = IntervalVector(x) - B2
  ctcPolar.contract(x_B2_r2_1[0],x_B2_r2_1[1], Interval(r2), Interval(theta_B2_r2_1))
  x_c |= (x_B2_r2_1 + B2)

  x_B3_r2_1 = IntervalVector(x) - B3
  ctcPolar.contract(x_B3_r2_1[0],x_B3_r2_1[1], Interval(r2), Interval(theta_B3_r2_1))
  x_c |= (x_B3_r2_1 + B3)

  
  x_B1_r4_1 = IntervalVector(x) - B1
  ctcPolar.contract(x_B1_r4_1[0],x_B1_r4_1[1], Interval(r4), Interval(theta_B1_r4_1))
  x_c |= (x_B1_r4_1 + B1)

  x_B2_r4_1 = IntervalVector(x) - B2
  ctcPolar.contract(x_B2_r4_1[0],x_B2_r4_1[1], Interval(r4), Interval(theta_B2_r4_1))
  x_c |= (x_B2_r4_1 + B2)

  x_B3_r4_1 = IntervalVector(x) - B3
  ctcPolar.contract(x_B3_r4_1[0],x_B3_r4_1[1], Interval(r4), Interval(theta_B3_r4_1))
  x_c |= (x_B3_r4_1 + B3)

  
  x_B1_r1_1 = IntervalVector(x) - B1
  ctcPolar.contract(x_B1_r1_1[0],x_B1_r1_1[1], Interval(r1), Interval(theta_B1_r1_1))
  x_c |= (x_B1_r1_1 + B1)

  x_B2_r1_1 = IntervalVector(x) - B2
  ctcPolar.contract(x_B2_r1_1[0],x_B2_r1_1[1], Interval(r1), Interval(theta_B2_r1_1))
  x_c |= (x_B2_r1_1 + B2)

  x_B3_r1_1 = IntervalVector(x) - B3
  ctcPolar.contract(x_B3_r1_1[0],x_B3_r1_1[1], Interval(r1), Interval(theta_B3_r1_1))
  x_c |= (x_B3_r1_1 + B3)


  x_A1_r3_1 = IntervalVector(x) - A1
  ctcPolar.contract(x_A1_r3_1[0],x_A1_r3_1[1], Interval(r3), Interval(theta_A1_r3_1))
  x_c |= (x_A1_r3_1 + A1)

  x_A2_r3_1 = IntervalVector(x) - A2
  ctcPolar.contract(x_A2_r3_1[0],x_A2_r3_1[1], Interval(r3), Interval(theta_A2_r3_1))
  x_c |= (x_A2_r3_1 + A2)

  x_A3_r3_1 = IntervalVector(x) - A3
  ctcPolar.contract(x_A3_r3_1[0],x_A3_r3_1[1], Interval(r3), Interval(theta_A3_r3_1))
  x_c |= (x_A3_r3_1 + A3)

  return x_c


def pave_contract(ctc, x, epsilon, X_plus):
  x = ctc(x)
  if not (x.is_empty()):
    if x.max_diam()<epsilon:
      X_plus.append(x)
    else:
      (x1,x2)=x.bisect_largest()
      pave_contract(ctc, x1, epsilon, X_plus)
      pave_contract(ctc, x2, epsilon, X_plus)

def pave_no_contract(ctc, x, epsilon, X_plus):
  xc = ctc(x) # the contraction is just used to assert that the constraint is satisfied, but the paving is done on the original box
  if not (xc.is_empty()):
    if x.max_diam()<epsilon:
      X_plus.append(x)
    else:
      (x1,x2)=x.bisect_largest()
      pave_no_contract(ctc, x1, epsilon, X_plus)
      pave_no_contract(ctc, x2, epsilon, X_plus)

# Points setup

O = Vector([0,0])

A1 = Vector ([0,5.35])
B1 = Vector ([0,3.9])

# Rotation Matrix of 2*PI/3
rot = Matrix([[np.cos(2*np.pi/3), np.sin(2*np.pi/3)],
              [-np.sin(2*np.pi/3), np.cos(2*np.pi/3)]])

A2 = rot*A1
B2 = rot*B1

A3 = rot*A2
B3 = rot*B2

# Distance between B1 and A2
d = (B1-A2).norm()

# Radiuses setup

R1, R2 = 9, 8.75

r1 = R1-B1[1]
r2 = R2-B1[1]
r3 = d - r1
r4 = 5.5
# r5 = r4 + 0.5
# r6 = R1-A1[1]

# Angles setup, notation is theta_center_radius_id

dtheta_0_R2= PI/2 + PI/40
theta_O_R2_1 = Interval([PI/2, PI/2 + dtheta_0_R2])
theta_O_R2_2 = Interval([PI + PI/6, PI+PI/6 + dtheta_0_R2])
theta_O_R2_3 = Interval([2*PI - PI/6, 2*PI - PI/6 + dtheta_0_R2])

theta_B1_r2_1 = Interval([-4.78*PI/12, PI/2])
theta_B2_r2_1 = theta_B1_r2_1 - 2*PI/3
theta_B3_r2_1 = theta_B2_r2_1 - 2*PI/3

theta_B1_r4_1 = Interval([-2.82*PI/12, 1.05*PI/4])
theta_B2_r4_1 = theta_B1_r4_1 - 2*PI/3
theta_B3_r4_1 = theta_B2_r4_1 - 2*PI/3

theta_B1_r1_1 = Interval([16.93*PI/12, 20.5*PI/12])
theta_B2_r1_1 = theta_B1_r1_1 - 2*PI/3
theta_B3_r1_1 = theta_B2_r1_1 - 2*PI/3

theta_A1_r3_1 = Interval([-8*PI/12, 14.75*PI/12])
theta_A2_r3_1 = theta_A1_r3_1 - 2*PI/3
theta_A3_r3_1 = theta_A2_r3_1 - 2*PI/3

# Setup for paving

x0 = IntervalVector([[-10,10],[-10,10]])

# Construction of the triskell

fig1 = Figure2D("Construction", GraphicOutput.VIBES)
fig1.set_window_properties([50,50],[500,500]) # set the window position and size
fig1.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

fig3 = Figure2D("Pave_without_contractors", GraphicOutput.VIBES | GraphicOutput.IPE)
fig3.set_window_properties([600,50],[500,500]) # set the window position and size
fig3.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1
fig3.draw_box(x0,StyleProperties.outside())

fig2 = Figure2D("Pave_with_contractors", GraphicOutput.VIBES | GraphicOutput.IPE)
fig2.set_window_properties([1150,50],[500,500]) # set the window position and size
fig2.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1
fig2.draw_box(x0,StyleProperties.outside())


## Support circles

fig1.draw_circle(O,R1,StyleProperties(".."))
fig1.draw_circle(O,R2,StyleProperties(".."))

fig1.draw_point(A1,StyleProperties(".."))
fig1.draw_point(B1,StyleProperties(".."))
fig1.draw_point(A2,StyleProperties(".."))
fig1.draw_point(B2,StyleProperties(".."))
fig1.draw_point(A3,StyleProperties(".."))
fig1.draw_point(B3,StyleProperties(".."))

fig1.draw_circle(B1,r1,StyleProperties(".."))
fig1.draw_circle(B2,r1,StyleProperties(".."))
fig1.draw_circle(B3,r1,StyleProperties(".."))
fig1.draw_circle(B1,r2,StyleProperties(".."))
fig1.draw_circle(B2,r2,StyleProperties(".."))
fig1.draw_circle(B3,r2,StyleProperties(".."))

fig1.draw_circle(A1,r3,StyleProperties(".."))
fig1.draw_circle(A2,r3,StyleProperties(".."))
fig1.draw_circle(A3,r3,StyleProperties(".."))

fig1.draw_circle(B1,r4,StyleProperties(".."))
fig1.draw_circle(B2,r4,StyleProperties(".."))
fig1.draw_circle(B3,r4,StyleProperties(".."))

## Support arcs

fig1.draw_pie(O,R2,theta_O_R2_1,StyleProperties(Color.red()))
fig1.draw_pie(O,R2,theta_O_R2_2,StyleProperties(Color.red()))
fig1.draw_pie(O,R2,theta_O_R2_3,StyleProperties(Color.red()))

fig1.draw_pie(B1,r2,theta_B1_r2_1,StyleProperties(Color.red()))
fig1.draw_pie(B2,r2,theta_B2_r2_1,StyleProperties(Color.red()))
fig1.draw_pie(B3,r2,theta_B3_r2_1,StyleProperties(Color.red()))

fig1.draw_pie(B1,r4,theta_B1_r4_1,StyleProperties(Color.red()))
fig1.draw_pie(B2,r4,theta_B2_r4_1,StyleProperties(Color.red()))
fig1.draw_pie(B3,r4,theta_B3_r4_1,StyleProperties(Color.red()))

fig1.draw_pie(B1,r1,theta_B1_r1_1,StyleProperties(Color.red()))
fig1.draw_pie(B2,r1,theta_B2_r1_1,StyleProperties(Color.red()))
fig1.draw_pie(B3,r1,theta_B3_r1_1,StyleProperties(Color.red()))

fig1.draw_pie(A1,r3,theta_A1_r3_1,StyleProperties(Color.red()))
fig1.draw_pie(A2,r3,theta_A2_r3_1,StyleProperties(Color.red()))
fig1.draw_pie(A3,r3,theta_A3_r3_1,StyleProperties(Color.red()))

# Paving

X_plus_no_contract = []
t0 = time.time()
pave_no_contract(CtcTryskell, IntervalVector(x0), 0.25, X_plus_no_contract)

area_no_contract = 0
for box in X_plus_no_contract:
  area_no_contract += box.volume()
  fig3.draw_box(box,StyleProperties.boundary())

print("Number of boxes in the paving without contractors: ", len(X_plus_no_contract))
print("Time taken for paving without contractors: ", time.time() - t0, " seconds")
print("Total area covered by the paving without contractors: ", area_no_contract)
print("Mean area of the boxes in the paving without contractors: ", area_no_contract/len(X_plus_no_contract))

X_plus_contract = []
t0 = time.time()
pave_contract(CtcTryskell, IntervalVector(x0), 0.25, X_plus_contract)

area_contract = 0
for box in X_plus_contract:
  area_contract += box.volume()
  fig2.draw_box(box,StyleProperties.boundary())

print("Number of boxes in the paving with contractors: ", len(X_plus_contract))
print("Time taken for paving with contractors: ", time.time() - t0, " seconds")
print("Total area covered by the paving with contractors: ", area_contract)
print("Mean area of the boxes in the paving with contractors: ", area_contract/len(X_plus_contract))
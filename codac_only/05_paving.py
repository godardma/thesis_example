from codac import *
import numpy as np

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
r4 = r1 + 0.27
# r5 = r4 + 0.5
# r6 = R1-A1[1]

# Angles setup, notation is theta_center_radius_id

dtheta_0_R2= PI/2 + PI/24
theta_O_R2_1 = Interval([PI/2, PI/2 + dtheta_0_R2])
theta_O_R2_2 = Interval([PI + PI/6, PI+PI/6 + dtheta_0_R2])
theta_O_R2_3 = Interval([2*PI - PI/6, 2*PI - PI/6 + dtheta_0_R2])

theta_B1_r2_1 = Interval([-4.75*PI/12, PI/2])
theta_B2_r2_1 = theta_B1_r2_1 - 2*PI/3
theta_B3_r2_1 = theta_B2_r2_1 - 2*PI/3

theta_B1_r4_1 = Interval([-2.95*PI/12, 1.14*PI/4])
theta_B2_r4_1 = theta_B1_r4_1 - 2*PI/3
theta_B3_r4_1 = theta_B2_r4_1 - 2*PI/3

theta_B1_r1_1 = Interval([16.95*PI/12, 20.5*PI/12])
theta_B2_r1_1 = theta_B1_r1_1 - 2*PI/3
theta_B3_r1_1 = theta_B2_r1_1 - 2*PI/3

theta_A1_r3_1 = Interval([-8*PI/12, 15*PI/12])
theta_A2_r3_1 = theta_A1_r3_1 - 2*PI/3
theta_A3_r3_1 = theta_A2_r3_1 - 2*PI/3

# Setup for paving

x0 = IntervalVector([[-10,10],[-10,10]])

# Construction of the triskell

fig1 = Figure2D("Construction", GraphicOutput.VIBES)
fig1.set_window_properties([50,50],[500,500]) # set the window position and size
fig1.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

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

fig1.draw_pie(O,R2,theta_O_R2_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(O,R2,theta_O_R2_2,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(O,R2,theta_O_R2_3,StyleProperties(Color.red(),"w:0.05"))

fig1.draw_pie(B1,r2,theta_B1_r2_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(B2,r2,theta_B2_r2_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(B3,r2,theta_B3_r2_1,StyleProperties(Color.red(),"w:0.05"))

fig1.draw_pie(B1,r4,theta_B1_r4_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(B2,r4,theta_B2_r4_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(B3,r4,theta_B3_r4_1,StyleProperties(Color.red(),"w:0.05"))

fig1.draw_pie(B1,r1,theta_B1_r1_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(B2,r1,theta_B2_r1_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(B3,r1,theta_B3_r1_1,StyleProperties(Color.red(),"w:0.05"))

fig1.draw_pie(A1,r3,theta_A1_r3_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(A2,r3,theta_A2_r3_1,StyleProperties(Color.red(),"w:0.05"))
fig1.draw_pie(A3,r3,theta_A3_r3_1,StyleProperties(Color.red(),"w:0.05"))

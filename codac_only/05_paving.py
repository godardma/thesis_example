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
r4 = r1 + 0.25
# r5 = r4 + 0.5
# r6 = R1-A1[1]

# Angles setup, notation is theta_center_radius_id

theta_O_R2_1 = Interval([PI/2,PI+PI/6])

# Setup fir paving

x0 = IntervalVector([[-10,10],[-10,10]])

# Construction of the triskell

fig1 = Figure2D("Construction", GraphicOutput.VIBES)
fig1.set_window_properties([50,50],[500,500]) # set the window position and size
fig1.set_axes(x0) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

## Support circles

fig1.draw_circle(O,R1)
fig1.draw_circle(O,R2)

fig1.draw_point(A1)
fig1.draw_point(B1)
fig1.draw_point(A2)
fig1.draw_point(B2)
fig1.draw_point(A3)
fig1.draw_point(B3)

fig1.draw_circle(B1,r1)
fig1.draw_circle(B2,r1)
fig1.draw_circle(B3,r1)
fig1.draw_circle(B1,r2)
fig1.draw_circle(B2,r2)
fig1.draw_circle(B3,r2)

fig1.draw_circle(A1,r3)
fig1.draw_circle(A2,r3)
fig1.draw_circle(A3,r3)

fig1.draw_circle(B2,r4)

## Support arcs

fig1.draw_pie(O,R2,theta_O_R2_1,StyleProperties(Color.red(),"w:0.1"))
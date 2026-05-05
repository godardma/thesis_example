from codac import *
import numpy as np

R1, R2 = 9, 8.75

A1 = Vector ([0,5.35])
B1 = Vector ([0,3.9])

# Rotation Matrix of 2*PI/3
rot = Matrix([[np.cos(2*np.pi/3), np.sin(2*np.pi/3)],
              [-np.sin(2*np.pi/3), np.cos(2*np.pi/3)]])

A2 = rot*A1
B2 = rot*B1

A3 = rot*A2
B3 = rot*B2

# distance between B1 and A2
d = (B1-A2).norm()

r1 = R1-B1[1]
r2 = R2-B1[1]
r3 = d - r1
r4 = r1 + 0.5
r5 = r4 + 0.5
r6 = R1-A1[1]

DefaultFigure.draw_circle([0,0],R1)
DefaultFigure.draw_circle([0,0],R2)

DefaultFigure.draw_point(A1)
DefaultFigure.draw_point(B1)
DefaultFigure.draw_point(A2)
DefaultFigure.draw_point(B2)
DefaultFigure.draw_point(A3)
DefaultFigure.draw_point(B3)

DefaultFigure.draw_circle(B1,r1)
DefaultFigure.draw_circle(B2,r1)
DefaultFigure.draw_circle(B3,r1)
DefaultFigure.draw_circle(B1,r2)
DefaultFigure.draw_circle(B2,r2)
DefaultFigure.draw_circle(B3,r2)

DefaultFigure.draw_circle(A1,r3)
DefaultFigure.draw_circle(A2,r3)
DefaultFigure.draw_circle(A3,r3)

DefaultFigure.draw_circle(A1,r6)
DefaultFigure.draw_circle(A2,r6)
DefaultFigure.draw_circle(A3,r6)

DefaultFigure.draw_circle(B2,r4)
DefaultFigure.draw_circle(B2,r5)

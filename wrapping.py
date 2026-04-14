from codac import *

fig = Figure2D("Wrapping rotation", GraphicOutput.VIBES | GraphicOutput.IPE)
fig.set_window_properties([50,50],[500,500]) # set the window position and size
fig.set_axes(axis(0,[-1.4,0.5]), axis(1,[-0.5,1.4])) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1

X0 = IntervalVector([[-0.1,0.1],[0.95,1.05]]) # initial set

# matrix rotation pi/4
mat1 = IntervalMatrix([[cos(PI/4), -sin(PI/4)],
            [sin(PI/4), cos(PI/4)]])

print(mat1)

# matrix rotation pi/2
mat2 = IntervalMatrix([[0,-1],[1,0]])


X1 = mat1*X0
X2_1 = mat1*X1
X2_2 = mat2*X0

test = X1 - X1.mid()
x_ub, y_ub = test[0].ub(), test[1].ub()
print(x_ub,y_ub)
Ab = sqrt(2)*Matrix([[y_ub, x_ub],[-y_ub, x_ub]])/2

z = IntervalVector([-sqrt(2)/2, sqrt(2)/2])
A = sqrt(2)*IntervalMatrix([[0.05, 0.1],[-0.05, 0.1]])/2

par = Parallelepiped(z.mid(),A.mid())
par2 = Parallelepiped(X2_2.mid(),Ab.mid())

# fig.draw_circle([0,0],1, Color.gray())
fig.draw_pie([0,0],Interval(1),Interval(PI/2,PI).inflate(0.3), Color.gray())
fig.draw_box(X0,[Color.black(),Color.light_blue()])
fig.draw_parallelepiped(par,StyleProperties([Color.black(),Color.light_blue()]))
fig.draw_parallelepiped(par2,StyleProperties([Color.dark_green()]))

fig.draw_box(X2_2,[Color.black(),Color.light_blue()])

fig.draw_box(X1,Color.red())
fig.draw_box(X2_1,Color.red())
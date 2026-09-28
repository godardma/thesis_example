# This is the example 04 of codac, explored_area (PEIBOS version)

from codac import *
import numpy as np

if __name__=="__main__":
  set_nb_threads(max_threads())

  # {psi0,Sigma} is a gnomonic atlas of the box [-1,1]^2
  X = VectorVar(1)
  psi0 = AnalyticFunction([X],[X[0],1])

  id = OctaSym([1,2])
  s = OctaSym([-2,1])

  Sigma = [id,s,s*s,s.invert()]

  # In this box, we consider that the x-axis is the width of the linear sensor
  # and the y-axis is the time
  y = VectorVar(2)
  L = 0.1*y[0]
  t = 1.2*y[1]

  # We construct the trajectory of the robot
  traj = vec(pow(t,3)-t,1-sqr(t))
  # We need its derivative to compute the orthogonal to the trajectory (for the sensor)
  dtraj = vec(3*sqr(t)-1, -2*t)
  dtraj_norm = sqrt(sqr(dtraj[0])+sqr(dtraj[1]))

  # The image of the box [-1,1]^2 by f is the swept area
  f = AnalyticFunction([y], [traj[0]-L*dtraj[1]/dtraj_norm, traj[1]+L*dtraj[0]/dtraj_norm])

  # For the SepImage, we need a contractor on the initial set (here a simple box)
  X0 = IntervalVector.constant(2,[-1,1])
  ctc_in = CtcWrapper(X0)

  # Separator on the area seen by a robot

  v_par = PEIBOS(f,psi0,Sigma,0.0625)

  sep = SepImage(f,psi0,Sigma,0.0625,ctc_in)

  # Visualizing the separator
  Y0 = IntervalVector([[-0.8,0.8],[-0.7,1.3]])

  fig_parallel = Figure2D("explored_area_parallel", GraphicOutput.VIBES|GraphicOutput.IPE) 
  fig_parallel.set_window_properties([50,50],[500,500])
  fig_parallel.set_axes(Y0)

  fig_pave = Figure2D("explored_area_pave", GraphicOutput.VIBES|GraphicOutput.IPE) 
  fig_pave.set_window_properties([600,50],[500,500])
  fig_pave.set_axes(Y0)

  for par in v_par :
    fig_parallel.draw_parallelepiped(par,StyleProperties.boundary())

  fig_pave.pave(Y0,sep,0.01)


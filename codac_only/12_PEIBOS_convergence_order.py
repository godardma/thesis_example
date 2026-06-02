import numpy as np
import matplotlib.pyplot as plt

epsilon = np.array([2**-2,2**-3,2**-4,2**-5,2**-6,2**-7,2**-8,2**-9])
vol_parallelepiped = np.array([12.366266808077224,1.2390858835282976,0.17580612575584198,0.03168810708178281,0.006629672829102125,0.0015098227331895205,0.00035984621146866976,8.78119045048233e-05])

# plot log10 of the volume of the parallelepiped against log10 of epsilon

plt.figure()
plt.plot(np.log(epsilon), np.log(vol_parallelepiped), marker='o')
plt.xlabel('log(epsilon)')
plt.ylabel('log(Volume of Parallelepiped)')
plt.show()
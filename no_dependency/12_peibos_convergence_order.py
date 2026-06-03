import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

epsilon = np.array([2**-2,2**-3,2**-4,2**-5,2**-6,2**-7,2**-8])
vol_box_natural = np.array([31.154720336490083,15.514955110738885,7.749346814248769,3.87364717412628,1.9366950018137488,0.9683314182174281,0.48416369849627566])
vol_box_centered = np.array([35.79089530176603,8.002908920156262,2.539069348532126,0.9936135040148106,0.4372319956647187,0.20479724380522918,0.0990726170923619])
vol_parallelepiped = np.array([12.366266808077224,1.2390858835282976,0.17580612575584198,0.03168810708178281,0.006629672829102125,0.0015098227331895205,0.00035984621146866976])




fig,ax = plt.subplots()

ax.plot(np.log(epsilon), np.log(vol_box_natural), marker='o', label='Natural box Enclosure')
ax.add_patch(Rectangle((np.log(epsilon[-1]), np.log(vol_box_natural[-1])), 0.5, 0.5, fill=False, edgecolor='blue', label='Natural box convergence order'))

ax.plot(np.log(epsilon), np.log(vol_parallelepiped), marker='o', label='Parallelepiped Enclosure')
ax.add_patch(Rectangle((np.log(epsilon[-1]), np.log(vol_parallelepiped[-1])), 0.5, 1, fill=False, edgecolor='orange', label='Parallelepiped convergence order'))

# ax.plot(np.log(epsilon), np.log(vol_box_centered), marker='o', label='Centered box Enclosure')
# ax.add_patch(Rectangle((np.log(epsilon[-1]), np.log(vol_box_centered[-1])), 0.5, 0.5, fill=False, edgecolor='green', label='Centered box convergence order'))

# ax.legend()
plt.show()
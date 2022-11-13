# Plot the bulk band structure, and compare it to the finite band structure.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider


a_x = 1   # Lattice constant
a_y = 1
m = 0.0  # Mass term (onsite energy difference)

kx = np.arange(-np.pi, np.pi, 0.005)
ky = np.arange(-np.pi, np.pi, 0.005)
kx, ky = np.meshgrid(kx, ky)

E = np.sqrt((m+np.cos(a_x*kx)+np.cos(a_y*ky))**2 + (np.sin(a_x*kx))**2 + (np.sin(a_y*ky))**2)
Em= -np.sqrt((m+np.cos(a_x*kx)+np.cos(a_y*ky))**2 + (np.sin(a_x*kx))**2 + (np.sin(a_y*ky))**2)


class HaldaneSurfaces:
  def __init__(self, ax):
    self.surf1 = ax.plot_surface(kx, ky, E, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
                       
    self.surf2 = ax.plot_surface(kx, ky, Em, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    self.ax = ax

  def update_surfaces(self, val):
      m = val
      self.surf1.remove()
      self.surf2.remove()
      E = np.sqrt((m+np.cos(a_x*kx)+np.cos(a_y*ky))**2 + (np.sin(a_x*kx))**2 + (np.sin(a_y*ky))**2)
      Em= -np.sqrt((m+np.cos(a_x*kx)+np.cos(a_y*ky))**2 + (np.sin(a_x*kx))**2 + (np.sin(a_y*ky))**2)
      
      self.surf1 = self.ax.plot_surface(kx, ky, E, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)                     
      self.surf2 = self.ax.plot_surface(kx, ky, Em, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(bottom=0.25)

haldane_surfaces = HaldaneSurfaces(ax)

ax.set_zlim(-2.0, 2.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')
fig.colorbar(haldane_surfaces.surf2, shrink=0.5, aspect=5)
slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slider = Slider(slider_ax, 'm', -2.0, 2.0, valinit=m)

def slider_changed(val):
  haldane_surfaces.update_surfaces(val)
  fig.canvas.draw_idle()
slider.on_changed(slider_changed)

plt.show()


# Notes: 
# Without symmetry protection
# = simplified graphene model (semi-metal) + mass term + interactions = break time-reversal and spatial inversion symmetry
# time symmetry: laws of physics stay the same when reversing time.
# spatial inversion symmetry: laws of physics stay the same when reflecting over an axis
#  
# Bernavik ch.1-2  TOPOLOGICAL INSULATORS AND TOPOLOGICAL SUPERCONDUCTORS
# Fritz Primer on Topological Insulators p26 eqn.2.26 derive to eqn.2.27 square lattice Haldane model
# s px py orbital model on a square lattice -> hopping between s-s, s-px, s-py, px-py for all sites (sigma, pi hoppings) 
# Haldane and then Fritz -> finite model hamilonian 
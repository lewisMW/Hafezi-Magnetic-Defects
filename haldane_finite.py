# Plot the bulk band structure, and compare it to the finite band structure.
from typing import Collection
from matplotlib.collections import PolyCollection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider


t = 1.0 # NN  hopping parameter
a = 1   # Lattice constant
phi = np.pi/2 # NNN phase
M = 0.2  # Mass term (onsite energy difference)
td= 0.05#M/(3*np.sqrt(3)) # NNN hopping

L = 19 # Length (equal in both dimensions)
N_sites = 2*L**2 # Need 2 sites for every N*N unit cells

np.set_printoptions(precision=1,edgeitems=30, linewidth=100000)

H = np.zeros((N_sites,N_sites), dtype=complex)
for i in range(N_sites):
  for j in range(N_sites):
    i_y = np.floor(i/(2*L)) # Site index in a1 lattice direction (row)
    i_x = i%(2*L)           # Site index in a2 lattice direction (row)
    j_y = np.floor(j/(2*L)) # Site index in a1 lattice direction (col)
    j_x = j%(2*L)           # Site index in a2 lattice direction (col)
    i_is_b = i%2       # is the row an 'b' site?
    j_is_b = j%2       # is the col an 'b' site?

    if i == j: # Same site (onsite self energy)
      H[i][j] = -M if i_is_b else M
    elif np.abs(i_x-j_x) == 1 and i_y-j_y == 0 or (i_y-j_y == 1 and i_x-j_x == -1 and not i_is_b) or (i_y-j_y == -1 and i_x-j_x == 1 and i_is_b): # Nearest-Neighbour Hopping
      H[i][j] = -t
    elif np.abs(i_x-j_x) == 2 and i_y-j_y == 0 or (i_y-j_y == 1 and i_x-j_x in (0,-2)) or (i_y-j_y == -1 and i_x-j_x in (0,2)): # Next-Nearest-Neighbour Hopping
      sign_ = 1
      if i_y-j_y == 0: # in the a1 direction
        sign_ = i_x-j_x == 2 and not i_is_b or i_x-j_x == -2 and i_is_b
      elif i_y-j_y == 1: # a2 direction, i is below j
        sign_ = i_x-j_x == -2 and not i_is_b or i_x-j_x == 0 and i_is_b
      elif i_y-j_y == -1: # a2 direction, i is above j
        sign_ = i_x-j_x == 0 and not i_is_b or i_x-j_x == 2 and i_is_b

      if not sign_: sign_ = -1

      phase = np.exp(sign_*phi*1j)
      H[i][j] = -td * phase
    else:
      H[i][j] = 0
    
  #print(H[i])

energies, eigenstates = np.linalg.eig(H)
    
indexed_energies = sorted(enumerate(energies), key=lambda x: np.real(x[1]))
# indexed_energies_ = list(map(lambda x: np.real(x[1]), indexed_energies))
# plt.scatter([0]*len(indexed_energies_),indexed_energies_)
# plt.show()

class HaldaneSurfaces:
  def __init__(self, ax):
    self.surf1 = None
    self.ax = ax
    self.update_surfaces(0)
  def update_surfaces(self, e):
      if self.surf1: self.surf1.remove()
      print(np.real(indexed_energies[e][1]))
      wave = eigenstates[:,indexed_energies[e][0]]
      Xs = []
      Ys = []
      Zs = []
      for i in range(2*L**2):
        i_y = np.floor(i/(2*L)) # Site index in a1 lattice direction (row)
        i_x = i%(2*L)           # Site index in a2 lattice direction (row) 
        i_is_b = i%2

        x = np.sqrt(3)*a * np.floor(i_x/2) + np.sqrt(3)/2*a*(i_y + i_is_b)  # v4 * ix + 
        y = 3/2*a*i_y + a/2 * i_is_b

        Xs.append(x)
        Ys.append(y)
        Zs.append(np.abs(wave[i])**2)
      self.ax.clear()
      self.surf1 = self.ax.plot_trisurf(np.array(Xs), np.array(Ys), np.array(Zs), cmap=cm.jet, linewidth=0)
      self.ax.set_box_aspect((np.ptp(Xs), np.ptp(Ys), np.ptp(Zs)*100*L))
      self.ax.scatter(Xs,Ys,c='k')


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(bottom=0.25)

wavefunction = HaldaneSurfaces(ax)

min_E = np.real(min(energies,key= lambda x: np.real(x)))
max_E = np.real(max(energies,key= lambda x: np.real(x)))
E_diff = max_E - min_E

slider_ax  = fig.add_axes([0.25, 0.05, 0.65, 0.03])
slider = Slider(slider_ax, "E'", min_E, max_E, valinit=min_E)

def slider_changed(val):
  e = int((val - min_E)/(max_E-min_E) * (len(energies)-1) )
  wavefunction.update_surfaces(e)
  #fig.canvas.draw_idle()
slider.on_changed(slider_changed)

plt.show()


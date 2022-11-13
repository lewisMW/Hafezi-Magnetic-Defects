# Plot the bulk band structure, and compare it to the finite band structure.
from typing import Collection
from matplotlib.collections import PolyCollection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider


t = 1.0  # NN  hopping parameter
a = 1    # Lattice constant
phi = np.pi/2 # NNN phase
M = 0.2  # Mass term (onsite energy difference)
td= 0.05 # M/(3*np.sqrt(3)) # NNN hopping

L = 4   # Length (equal in both dimensions)
N_sites = 2*L**2 # Need 2 sites for every N*N unit cells

np.set_printoptions(precision=1,edgeitems=30, linewidth=100000)

class HaldaneSurfaces:
  def __init__(self, ax, ax2):
    self.ax = ax
    self.ax2 = ax2
    self.ax2_xlim = None
    self.ax2_ylim = None
    self.indexed_energies = None
    self.eigenstates = None
    self.e = None
    self.c1_i =None
    self.c1_j =None
    self.c1_ab =None
    self.t1_c =None
    self.c2_i =None
    self.c2_j =None
    self.c2_ab =None
    self.t2_c =None
    self.solve(0,0,0,0,0,0,0,0)
    self.update_surfaces(0)
    
  def update_surfaces(self, e=None):
      if e!=None: self.e = e
      print(np.real(self.indexed_energies[self.e][1]))
      wave = self.eigenstates[:,self.indexed_energies[self.e][0]]
      Xs = []
      Ys = []
      Zs_1 = []
      Zs_2 = []
      for i in range(2*L**2):
        i_y = np.floor(i/(2*L)) # Site index in a1 lattice direction (row)
        i_x = i%(2*L)           # Site index in a2 lattice direction (row) 
        i_is_b = i%2

        x = np.sqrt(3)*a * np.floor(i_x/2) + np.sqrt(3)/2*a*(i_y + i_is_b)  # v4 * ix + 
        y = 3/2*a*i_y + a/2 * i_is_b

        Xs.append(x)
        Ys.append(y)
        Zs_1.append(np.abs(wave[i])**2)
        Zs_2.append(np.abs(wave[i + 2*L**2])**2)
      self.ax.clear()
      #self.ax2.clear()
      #self.surf1 = self.ax.plot_trisurf(np.array(Xs), np.array(Ys), np.array(Zs_1), cmap=cm.jet, linewidth=0)
      #self.surf2 = self.ax2.plot_trisurf(np.array(Xs), np.array(Ys), np.array(Zs_2), cmap=cm.jet, linewidth=0)
      #self.ax.set_box_aspect((np.ptp(Xs), np.ptp(Ys), np.ptp(Zs_1)*100*L))
      self.ax.scatter(Xs,Ys,s=80, c=Zs_1,cmap=cm.jet)
      self.ax.set_xlim((-1,3*L))
      self.ax.set_ylim((-1,1.5*L))
      #self.ax2.set_box_aspect((np.ptp(Xs), np.ptp(Ys), np.ptp(Zs_1)*100*L))
      #self.ax2.scatter(Xs,Ys,c='k')
      energies = [np.real(energy[1]) for energy in self.indexed_energies]
      #self.ax2_xlim = self.ax2.get_xlim()
      #self.ax2_ylim = self.ax2.get_ylim()
      self.ax2.clear()
      self.ax2.scatter(Xs,Ys,s=80, c=Zs_2,cmap=cm.jet)
      self.ax2.set_xlim((-1,3*L))
      self.ax2.set_ylim((-1,1.5*L))
      #self.ax2.scatter(np.zeros((1,len(energies))),energies,marker='x',c="red",alpha=0.5,linewidths=0.5)
      #self.ax2.grid()
      #self.ax2.set_xlim(self.ax2_xlim)
      #self.ax2.set_ylim(self.ax2_ylim)

  def solve(self, c1_i=None, c1_j=None, c1_ab=None, t1_c=None, c2_i=None, c2_j=None, c2_ab=None, t2_c=None):
    if c1_i!=None:  self.c1_i =c1_i
    if c1_j!=None:  self.c1_j =c1_j
    if c1_ab!=None: self.c1_ab=c1_ab
    if t1_c!=None:  self.t1_c =t1_c
    if c2_i!=None:  self.c2_i =c2_i
    if c2_j!=None:  self.c2_j =c2_j
    if c2_ab!=None: self.c2_ab=c2_ab
    if t2_c!=None:  self.t2_c =t2_c
    # coupling_site is (0,0,0) to (L-1,L-1,1)
    coupling_index1 = self.c1_i * L*2 + self.c1_j* 2 + self.c1_ab
    coupling_index2 = self.c2_i * L*2 + self.c2_j* 2 + self.c2_ab
    H_c  = np.zeros((N_sites,N_sites), dtype=complex)
    for i in range(N_sites):
      for j in range(N_sites):
        i_y = np.floor(i/(2*L)) # Site index in a1 lattice direction (row)
        i_x = i%(2*L)           # Site index in a2 lattice direction (row)
        j_y = np.floor(j/(2*L)) # Site index in a1 lattice direction (col)
        j_x = j%(2*L)           # Site index in a2 lattice direction (col)
        i_is_b = i%2       # is the row an 'b' site?
        j_is_b = j%2       # is the col an 'b' site?

        if i == j: # Same site (onsite self energy)
          H_c[i][j] = -M if i_is_b else M
        elif np.abs(i_x-j_x) == 1 and i_y-j_y == 0 or (i_y-j_y == 1 and i_x-j_x == -1 and not i_is_b) or (i_y-j_y == -1 and i_x-j_x == 1 and i_is_b): # Nearest-Neighbour Hopping
          H_c[i][j] = -t
        elif np.abs(i_x-j_x) == 2 and i_y-j_y == 0 or (i_y-j_y == 1 and i_x-j_x in (0,-2)) or (i_y-j_y == -1 and i_x-j_x in (0,2)): # Next-Nearest-Neighbour Hopping
          sign_ = 1
          if i_y-j_y == 0: # in the a1 direction
            sign_ = i_x-j_x == 2 and not i_is_b or i_x-j_x == -2 and i_is_b
          elif i_y-j_y == 1: # a2 direction, i is below j
            sign_ = i_x-j_x == -2 and not i_is_b or i_x-j_x == 0 and i_is_b
          elif i_y-j_y == -1: # a2 direction, i is above j
            sign_ = i_x-j_x == 0 and not i_is_b or i_x-j_x == 2 and i_is_b

          if not sign_: sign_ = -1

          phase      =  np.exp(sign_*phi*1j)
          H_c[i][j]  = -td * phase
        else:
          H_c[i][j]  = 0
        
      #print(H[i])

    ### Coupling:

    C_  = np.zeros((N_sites,N_sites), dtype=complex)
    C_[coupling_index1][coupling_index1] = self.t1_c
    C_[coupling_index2][coupling_index2] = self.t2_c

    H = np.block([ [H_c , C_  ], [C_ , H_c.conj() ] ])

    np.savetxt('mat.csv' , H, delimiter=',', fmt='%s')

    energies, self.eigenstates = np.linalg.eig(H)
        
    self.indexed_energies = sorted(enumerate(energies), key=lambda x: np.real(x[1]))
    # indexed_energies_ = list(map(lambda x: np.real(x[1]), indexed_energies))
    # plt.scatter([0]*len(indexed_energies_),indexed_energies_)
    # plt.show()


fig, ax = plt.subplots()#subplot_kw={"projection": "3d"})
ax = fig.add_subplot(121)#, projection='3d')
ax2 = fig.add_subplot(122)#, projection='3d')
fig.subplots_adjust(bottom=0.4)

wavefunction = HaldaneSurfaces(ax, ax2)

slider_ax0  = fig.add_axes([0.15, 0.16, 0.70, 0.02])
slider0 = Slider(slider_ax0, "E'", int(N_sites*0.75),int(N_sites*1.25), valinit=wavefunction.e)    #2*N_sites-1, valinit=N_sites-1)

slider_ax1  = fig.add_axes([0.15, 0.14, 0.70, 0.02])
slider1 = Slider(slider_ax1, "c1_i", 0, L-1, valinit=wavefunction.c1_i)

slider_ax2  = fig.add_axes([0.15, 0.12, 0.70, 0.02])
slider2 = Slider(slider_ax2, "c1_j", 0, L-1, valinit=wavefunction.c1_j)

slider_ax3  = fig.add_axes([0.15, 0.10, 0.70, 0.02])
slider3 = Slider(slider_ax3, "c1_ab", 0, 1, valinit=wavefunction.c1_ab)

slider_ax4  = fig.add_axes([0.15, 0.08, 0.70, 0.02])
slider4 = Slider(slider_ax4, "t1_c", -4.0, 4.0, valinit=wavefunction.t1_c)

slider_ax5  = fig.add_axes([0.15, 0.06, 0.70, 0.02])
slider5 = Slider(slider_ax5, "c2_i", 0, L-1, valinit=wavefunction.c2_i)

slider_ax6  = fig.add_axes([0.15, 0.04, 0.70, 0.02])
slider6 = Slider(slider_ax6, "c2_j", 0, L-1, valinit=wavefunction.c2_j)

slider_ax7  = fig.add_axes([0.15, 0.02, 0.70, 0.02])
slider7 = Slider(slider_ax7, "c2_ab", 0, 1, valinit=wavefunction.c2_ab)

slider_ax8  = fig.add_axes([0.15, 0.00, 0.70, 0.02])
slider8 = Slider(slider_ax8, "t2_c", -4.0, 4.0, valinit=wavefunction.t2_c)

def slider_changed0(val):
  wavefunction.update_surfaces(int(val))
slider0.on_changed(slider_changed0)

def slider_changed1(val):
  wavefunction.solve(c1_i=int(val))
  wavefunction.update_surfaces()
slider1.on_changed(slider_changed1)

def slider_changed2(val):
  wavefunction.solve(c1_j=int(val))
  wavefunction.update_surfaces()
slider2.on_changed(slider_changed2)

def slider_changed3(val):
  wavefunction.solve(c1_ab=int(val))
  wavefunction.update_surfaces()
slider3.on_changed(slider_changed3)

def slider_changed4(val):
  wavefunction.solve(t1_c=val)
  wavefunction.update_surfaces()
slider4.on_changed(slider_changed4)

def slider_changed5(val):
  wavefunction.solve(c2_i=int(val))
  wavefunction.update_surfaces()
slider5.on_changed(slider_changed5)

def slider_changed6(val):
  wavefunction.solve(c2_j=int(val))
  wavefunction.update_surfaces()
slider6.on_changed(slider_changed6)

def slider_changed7(val):
  wavefunction.solve(c2_ab=int(val))
  wavefunction.update_surfaces()
slider7.on_changed(slider_changed7)

def slider_changed8(val):
  wavefunction.solve(t2_c=val)
  wavefunction.update_surfaces()
slider8.on_changed(slider_changed8)

plt.show()


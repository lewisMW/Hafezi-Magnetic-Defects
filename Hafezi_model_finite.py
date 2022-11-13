# Plot the bulk band structure, and compare it to the finite band structure.
from typing import Collection
from matplotlib.collections import PolyCollection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider, Button
import matplotlib.path as mpath

import os

def GenerateTrivialHamiltonian(Nx,Ny,d = 1.0,t = 0.1,t_= 0.5,K = 8.0):
    H = np.zeros([2*Nx*Ny, 2*Nx*Ny])
    for i in range(int(2*Nx*Ny)):
        for j in range(int(2*Nx*Ny)):
            if np.abs(i-j) == 1:
                if i>j and i%2==1 or j>i and j%2==1: H[i][j] = -t_*K
            elif np.abs(i-j) == 2 and np.floor(i/(2*Nx)) == np.floor(j/(2*Nx)):
                if i>j and i%2==0 or j>i and j%2==0: H[i][j] = -t*K
            elif i%2==0 and (j==i+2*Nx or i==j+2*Nx):
                H[i][j] = -t*K
    return H

alpha = 1.0/3.0
k = 1.0

a = 1.0

L = 21   # Length (equal in both dimensions)
N_sites = L**2 # Need 1 site for every N*N unit cells

T_rows=2

np.set_printoptions(precision=4,edgeitems=30, linewidth=100000)

circle0 = mpath.Path.circle(radius=1.0)
circle1 = mpath.Path.circle(radius=0.3)
ring_marker = mpath.Path(
    vertices=np.concatenate([circle0.vertices, circle1.vertices[::-1, ...]]),
    codes=np.concatenate([circle0.codes, circle1.codes]))

class HafeziSurface:
  def __init__(self, ax_spin1, ax_spin2, ax_energies, ax_couplings):
    self.folder_n = 0
    self.test_n = 0
    self.high_score = 0
    self.best_C = None#np.diag(np.loadtxt("C:\Users\lewis\OneDrive - The University of Sydney (Students)\Thesis_1_1\Work\C_climb\C_0.002518.csv", delimiter=','))
    self.surf1 = None
    self.surf2 = None
    self.ax_spin1 = ax_spin1
    self.ax_spin2 = ax_spin2
    self.ax_energies = ax_energies
    self.ax_couplings = ax_couplings
    self.ax_energies_xlim = None
    self.ax_energies_ylim = None
    self.indexed_energies = None
    self.eigenstates = None
    self.energies_N = 0
    self.e = None
    self.c1_i =None
    self.c1_j =None
    self.t1_c =None
    self.c2_i =None
    self.c2_j =None
    self.t2_c =None
    self.t =0.1
    self.t_ =0.5
    self.K = 8.0
    self.T = 1.0
    self.C = None
    self.solve(4,6,4,4,14,4,0,0,8,0)
    #self.solve(4,6,4,4,14,4,0.1,0.5,8,5)
    self.update_surfaces(671) # E
    ax_energies.set_xlim((-0.1,0.1))
    ax_energies.set_ylim((0.6,2))#(np.real(self.indexed_energies[0][1]),np.real(self.indexed_energies[-1][1])))
    
  def update_surfaces(self, e=None, measure=False):
      if e!=None: self.e = e
      this_energy = np.real(self.indexed_energies[self.e][1])
      print(this_energy)
      wave = self.eigenstates[:,self.indexed_energies[self.e][0]]
      self.energies_N = len(wave)
      Xs = []
      Ys = []
      Zs = []
      couple = np.diag(self.C)
      couplings = []
      pos = 0
      for i in range(int(len(wave)/2)):
        if i < 2*T_rows*L and i%2==1: continue # Skip B sites in the trivial insulator
        i_y = np.floor(pos/L)
        i_x = pos%(L)          

        x = a*i_x
        y = a*i_y

        Xs.append(x)
        Ys.append(y)
        Zs.append(np.abs(wave[i])**2)
        couplings.append(couple[i])
        pos += 1
      for i in range(int(len(wave)/2), len(wave)):
        if i < int(len(wave)/2)+2*T_rows*L and (i-int(len(wave)/2))%2==1: continue # Skip B sites in the trivial insulator
        Zs.append(np.abs(wave[i])**2)
      Zs = np.asarray(Zs)
      Zs_norm = np.sqrt((Zs).sum())
      Zs = Zs/Zs_norm if Zs_norm != 0 else Zs*0
      Zs_1 = Zs[:int(len(Zs)/2) ]
      Zs_2 = Zs[ int(len(Zs)/2):]

      normalise = colors.Normalize(vmin=0, vmax=np.max(Zs))
      self.ax_spin1.clear()
      self.ax_spin1.scatter(Xs,Ys,s=80, c=Zs_1,norm=normalise,cmap=cm.jet, marker=ring_marker,edgecolors='black',linewidth=0.2)
      self.ax_spin2.clear()
      self.ax_spin2.scatter(Xs,Ys,s=80, c=Zs_2,norm=normalise,cmap=cm.jet, marker=ring_marker,edgecolors='black',linewidth=0.2)

      energies = [np.real(energy[1]) for energy in self.indexed_energies]
      self.ax_energies_xlim = self.ax_energies.get_xlim()
      self.ax_energies_ylim = self.ax_energies.get_ylim()
      self.ax_energies.clear()
      self.ax_energies.scatter(np.zeros((1,len(energies))),energies,marker='x',c="red",alpha=0.5,linewidths=0.5)
      self.ax_energies.scatter(0,np.real(this_energy),marker='x',c="blue",alpha=1,linewidths=0.5)
      self.ax_energies.grid()
      self.ax_energies.set_xlim(self.ax_energies_xlim)
      self.ax_energies.set_ylim(self.ax_energies_ylim)

      self.ax_couplings.clear()
      self.ax_couplings.scatter(Xs,Ys,s=80, c=couplings,cmap=cm.jet, marker=ring_marker,edgecolors='black',linewidth=0.2)

      if measure:
        M1 = np.reshape(Zs_1, (-1,L))
        M1 = M1[T_rows:, 0:]
        edges_score = np.min(np.concatenate((M1[0,0:L],M1[L-1,0:L],M1[1:L-1,0],M1[1:L-1,L-1])))
        print("Minimum edge site:",edges_score)
        #print("Score:",edges_score,"compared to high score",self.high_score)
        #if edges_score > self.high_score:
        #  self.high_score = edges_score
        #  self.best_C = self.C
        #  np.savetxt('C_climb/%d/C_%f.csv'%(self.folder_n,edges_score), np.real(np.diag(self.C)), delimiter=',', fmt='%s')
        #  plt.savefig('C_climb/%d/state_%f.png'%(self.folder_n,edges_score))
        
        # if self.test_n % 1000 == 0:
        #   self.best_C = None
        #   self.C = None
        #   self.high_score = 0
        # if self.test_n % 5000 == 0:
        #   self.folder_n += 1
        #   while os.path.isdir("C_climb/%d" % self.folder_n): self.folder_n += 1
        #   os.mkdir("C_climb/%d"%self.folder_n)

        #   edge_pick = np.random.randint(4)
        #   if edge_pick == 0:
        #     self.c1_i = 0
        #     self.c1_j = np.random.randint(L)
        #   if edge_pick == 1:
        #     self.c1_i = np.random.randint(1,L-1)
        #     self.c1_j = 0
        #   if edge_pick == 2:
        #     self.c1_i = L-1
        #     self.c1_j = np.random.randint(L)
        #   if edge_pick == 3:
        #     self.c1_i = np.random.randint(1,L-1)
        #     self.c1_j = L-1
        #   edge_pick = np.random.randint(4)
        #   if edge_pick == 0:
        #     self.c2_i = 0
        #     self.c2_j = np.random.randint(L)
        #   if edge_pick == 1:
        #     self.c2_i = np.random.randint(1,L-1)
        #     self.c2_j = 0
        #   if edge_pick == 2:
        #     self.c2_i = L-1
        #     self.c2_j = np.random.randint(L)
        #   if edge_pick == 3:
        #     self.c2_i = np.random.randint(1,L-1)
        #     self.c2_j = L-1

        # self.test_n += 1  

  def solve(self, c1_i=None, c1_j=None, t1_c=None, c2_i=None, c2_j=None, t2_c=None, t=None,t_=None,K=None, T=None,use_best_C=False):
    if c1_i!=None:  self.c1_i =c1_i
    if c1_j!=None:  self.c1_j =c1_j
    if t1_c!=None:  self.t1_c =t1_c
    if c2_i!=None:  self.c2_i =c2_i
    if c2_j!=None:  self.c2_j =c2_j
    if t2_c!=None:  self.t2_c =t2_c
    if t!=None:  self.t = t
    if t_!=None:  self.t_ = t_
    if K!=None:  self.K = K
    if T!=None:   self.T = T
    # coupling_site is (0,0,0) to (L-1,L-1,1)
    coupling_index1 = self.c1_i * L + self.c1_j
    coupling_index2 = self.c2_i * L + self.c2_j

    # Generate Hamiltonian:
    H_c  = np.zeros((N_sites,N_sites), dtype=complex)
    for i in range(N_sites):
      for j in range(N_sites):
        if abs(i-j) == 1 and int(i/L) == int(j/L): # just-off-diagonal terms. Check these are not edges (adjacent sites in the same y -> center block-diagonal)
            y = int(i/L)
            phase =  2*np.pi*alpha*1j*y
            if i > j: # negtive phase
              phase = -phase
            H_c[i][j] = -k*np.exp(phase)  
        elif abs(i-j) == L: # vertically adjacent sites
          H_c[i][j] = -k
        else:
          H_c[i][j] = 0
    H_t = GenerateTrivialHamiltonian(L,T_rows,1,self.t,self.t_,self.K)
    T_ = np.zeros((2*T_rows*L, N_sites))
    for i in range(L):
      T_[L*T_rows+2*i,i] = -self.T
    H_c = np.block([[H_t, T_], [T_.conj().transpose(),H_c]])
    #np.savetxt('mat.csv', H_c, delimiter=',', fmt='%s')

    ### Coupling:
    #if use_best_C:
    #if type(self.best_C) != type(None): self.C = self.best_C.copy()
    #else:
    #if type(self.C) == type(None):
    self.C  = np.zeros(np.shape(H_c), dtype=complex)
    # for i in range(np.random.randint(10)):
    #   rand_site = np.random.randint(N_sites)
    #   self.C[rand_site][rand_site] = np.random.random()*16-8# * np.exp(np.random.random()*2*np.pi)
    # Completely random:
    #for i in range(np.shape(self.C)[0]):#2*T_rows*L,2*T_rows*L+L):#np.shape(self.C)[0]):
    #  self.C[i][i] = 0.005*np.random.random()*np.exp(np.random.random()*2*np.pi)
    self.C[coupling_index1][coupling_index1] = self.t1_c
    self.C[coupling_index2][coupling_index2] = self.t2_c

    H = np.block([ [H_c , self.C  ], [self.C.conj().transpose() , H_c.conj() ] ])

    energies, self.eigenstates = np.linalg.eig(H)
        
    self.indexed_energies = sorted(enumerate(energies), key=lambda x: np.real(x[1]))






fig, ax = plt.subplots()#subplot_kw={"projection": "3d"})
ax_spin1 = fig.add_subplot(141)#, projection='3d')
ax_spin2 = fig.add_subplot(142)
ax_energies = fig.add_subplot(143)
ax_couplings = fig.add_subplot(144)
fig.subplots_adjust(bottom=0.4)

wavefunction = HafeziSurface(ax_spin1,ax_spin2,ax_energies,ax_couplings)

gbutton_ax = fig.add_axes([0.15, 0.24, 0.70, 0.02])
gbutton    = Button(gbutton_ax, "Generate")

sbutton_ax = fig.add_axes([0.15, 0.22, 0.70, 0.02])
sbutton    = Button(sbutton_ax, "Save Best")

E_indexs = np.linspace(0,wavefunction.energies_N, wavefunction.energies_N)
slider_ax0  = fig.add_axes([0.15, 0.20, 0.70, 0.02])
slider0 = Slider(slider_ax0, "E'", 0, wavefunction.energies_N, valinit=wavefunction.e, valstep=E_indexs)#0, 2*N_sites-1, valinit=wavefunction.e, valstep=E_indexs)

slider_ax1  = fig.add_axes([0.15, 0.18, 0.70, 0.02])
slider1 = Slider(slider_ax1, "c1_i", 0, L-1, valinit=wavefunction.c1_i)

slider_ax2  = fig.add_axes([0.15, 0.16, 0.70, 0.02])
slider2 = Slider(slider_ax2, "c1_j", 0, L-1, valinit=wavefunction.c1_j)

slider_ax3  = fig.add_axes([0.15, 0.14, 0.70, 0.02])
slider3 = Slider(slider_ax3, "t1_c", -4.0, 20.0, valinit=wavefunction.t1_c)

slider_ax4  = fig.add_axes([0.15, 0.12, 0.70, 0.02])
slider4 = Slider(slider_ax4, "c2_i", 0, L-1, valinit=wavefunction.c2_i)

slider_ax5  = fig.add_axes([0.15, 0.10, 0.70, 0.02])
slider5 = Slider(slider_ax5, "c2_j", 0, L-1, valinit=wavefunction.c2_j)

slider_ax6  = fig.add_axes([0.15, 0.08, 0.70, 0.02])
slider6 = Slider(slider_ax6, "t2_c", -4.0, 20.0, valinit=wavefunction.t2_c)

slider_ax7  = fig.add_axes([0.15, 0.06, 0.70, 0.02])
slider7 = Slider(slider_ax7, "t", -2, 2, valinit=wavefunction.t)

slider_ax8  = fig.add_axes([0.15, 0.04, 0.70, 0.02])
slider8 = Slider(slider_ax8, "t_", -2, 2, valinit=wavefunction.t_)

slider_ax9  = fig.add_axes([0.15, 0.02, 0.70, 0.02])
slider9 = Slider(slider_ax9, "K", -16, 16, valinit=wavefunction.K)

slider_ax10  = fig.add_axes([0.15, 0.00, 0.70, 0.02])
slider10 = Slider(slider_ax10, "T", -5, 5, valinit=wavefunction.T)

def gbutton_clicked(self):
  #while True:
    wavefunction.solve()
    wavefunction.update_surfaces(measure=True)
    
  #t2_c=wavefunction.t2_c+0.2)
  #print("Changed t2_c to %f" % wavefunction.t2_c)
  #wavefunction.solve(c1_i=wavefunction.c1_i+1)
  #print("Changed c1_i to %f" % wavefunction.c1_i)
gbutton.on_clicked(gbutton_clicked)

def sbutton_clicked(self):
  wavefunction.solve(use_best_C=True)
  wavefunction.update_surfaces()

sbutton.on_clicked(sbutton_clicked)

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
  wavefunction.solve(t1_c=val)
  wavefunction.update_surfaces()
slider3.on_changed(slider_changed3)

def slider_changed4(val):
  wavefunction.solve(c2_i=int(val))
  wavefunction.update_surfaces()
slider4.on_changed(slider_changed4)

def slider_changed5(val):
  wavefunction.solve(c2_j=int(val))
  wavefunction.update_surfaces()
slider5.on_changed(slider_changed5)

def slider_changed6(val):
  wavefunction.solve(t2_c=val)
  wavefunction.update_surfaces()
slider6.on_changed(slider_changed6)

def slider_changed7(val):
  wavefunction.solve(t=val)
  wavefunction.update_surfaces()
slider7.on_changed(slider_changed7)

def slider_changed8(val):
  wavefunction.solve(t_=val)
  wavefunction.update_surfaces()
slider8.on_changed(slider_changed8)

def slider_changed9(val):
  wavefunction.solve(K=val)
  wavefunction.update_surfaces()
slider9.on_changed(slider_changed9)

def slider_changed10(val):
  wavefunction.solve(T=val)
  wavefunction.update_surfaces()
slider10.on_changed(slider_changed10)

plt.show()

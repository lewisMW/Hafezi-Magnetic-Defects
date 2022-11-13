# Plot the bulk band structure, and compare it to the finite band structure.
import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider

np.set_printoptions(precision=1,edgeitems=30, linewidth=100000)

kx_ = np.arange(-np.pi, np.pi*1.01, np.pi/100)
N = 21 # sites in the y direction
d = 1.0
t = 0.1
t_= 0.5
K = 8.0

def GenerateHamiltonian(N = 21,d = 1.0,t = 0.1,t_= 0.5,K = 8.0):
    H = np.zeros([2*N*N, 2*N*N])
    for i in range(int(2*N*N)):
        for j in range(int(2*N*N)):
            if np.abs(i-j) == 1:
                if i>j and i%2==1 or j>i and j%2==1: H[i][j] = -t_*K
            elif np.abs(i-j) == 2:
                if i>j and i%2==0 or j>i and j%2==0: H[i][j] = -t*K
            elif i%2==0 and (j==i+2*N or i==j+2*N):
                H[i][j] = -t*K
    return H

energies, eigenstates = np.linalg.eig(GenerateHamiltonian())
plot.scatter(np.zeros((2*N*N)), np.real(energies))
plot.show()
input("hi")

def CalcSpectra():
  plotData = np.zeros([len(kx_),2*N])
  for kx_i, kx in enumerate(kx_):
    H = np.zeros([2*N, 2*N])
    for i in range(int(2*N)):
      for j in range(int(2*N)):
        if np.abs(i-j) == 1:
          if i>j and i%2==1 or j>i and j%2==1: H[i][j] = -t_*K
        if np.abs(i-j) == 2:
          if i>j and i%2==0 or j>i and j%2==0: H[i][j] = -t*K
        elif i == j:
          if i%2==0:
              H[i][j] = -2*t*np.cos(kx*d)*K
    energies, eigenstates = np.linalg.eig(H)
    plotData[kx_i] = np.sort(energies)
  return plotData


fig = plot.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3)
lines = ax.plot(kx_,CalcSpectra())
plot.title("Trivial Insulator Model Energies")
plot.xlabel("$k_x$")
plot.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],["$-\\pi$","$-\\frac{\\pi}{2}$","$0$","$\\frac{\\pi}{2}$","$\\pi$"])
plot.ylabel("Energy")
slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slider = Slider(slider_ax, '$m$', -2.5, 2.5, valinit=1)
slider_td_ax  = fig.add_axes([0.25, 0.05, 0.65, 0.03])
slider_td = Slider(slider_td_ax, '$t_x$', 0.0, 0.5, valinit=1)



def slider_changed(val):
    lines_data = CalcSpectra()
    for i,l in enumerate(lines):
      l.set_ydata(lines_data[:,i])
    fig.canvas.draw_idle()
slider.on_changed(slider_changed)

def slider_td_changed(val):
    lines_data = CalcSpectra()
    for i,l in enumerate(lines):
      l.set_ydata(lines_data[:,i])
    fig.canvas.draw_idle()
slider_td.on_changed(slider_td_changed)

plot.show()

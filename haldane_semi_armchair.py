# Plot the bulk band structure, and compare it to the finite band structure.
import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider

np.set_printoptions(precision=1,edgeitems=30, linewidth=100000)

a_x = 1   # Lattice constant
m = 0  # Mass term (onsite energy difference)
t = 1.0 # NN  hopping parameter
phi = np.pi/2 # NNN phase
td= 0.0#M/(3*np.sqrt(3)) # NNN hopping
kx_ = np.arange(-np.pi, np.pi*1.01, np.pi/100)
N = 10 # sites in the y direction

def CalcSpectra(m, td):
  plotData = np.zeros([len(kx_),2*N])
  for kx_i, kx in enumerate(kx_):
    H = np.zeros([2*N, 2*N])
    for i in range(int(N)):
      for j in range(int(N)):
        if i == j:
          H[i*2][j*2] = m - 2*td* np.cos(np.sqrt(3)*a_x*kx-phi)
          H[i*2+1][j*2] = -2*t*np.cos(np.sqrt(3)/2*a_x*kx)
          H[i*2][j*2+1] = -2*t*np.cos(np.sqrt(3)/2*a_x*kx)
          H[i*2+1][j*2+1] = -m - 2*td* np.cos(np.sqrt(3)*a_x*kx+phi)
        if i == j+1:
          H[i*2][j*2] = -2*td*np.cos(np.sqrt(3)/2*a_x*kx+phi)
          H[i*2+1][j*2] = -t
          H[i*2][j*2+1] = 0
          H[i*2+1][j*2+1] = -2*td*np.cos(np.sqrt(3)/2*a_x*kx-phi)
        if i+1 == j:
          H[i*2][j*2] = -2*td*np.cos(np.sqrt(3)/2*a_x*kx+phi)
          H[i*2+1][j*2] = 0
          H[i*2][j*2+1] = -t
          H[i*2+1][j*2+1] = -2*td*np.cos(np.sqrt(3)/2*a_x*kx-phi)
    energies, eigenstates = np.linalg.eig(H)
    plotData[kx_i] = np.sort(energies)
  return plotData


fig = plot.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3)
lines = ax.plot(kx_,CalcSpectra(m,td))
plot.title("Haldane Semi-Infinite Model Energies")
plot.xlabel("$k_x$")
plot.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],["$-\\pi$","$-\\frac{\\pi}{2}$","$0$","$\\frac{\\pi}{2}$","$\\pi$"])
plot.ylabel("Energy")
slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slider = Slider(slider_ax, '$m$', -2.5, 2.5, valinit=m)
slider_td_ax  = fig.add_axes([0.25, 0.05, 0.65, 0.03])
slider_td = Slider(slider_td_ax, '$t_x$', 0.0, 0.5, valinit=td)



def slider_changed(val):
    lines_data = CalcSpectra(m=val, td=slider_td.val)
    for i,l in enumerate(lines):
      l.set_ydata(lines_data[:,i])
    fig.canvas.draw_idle()
slider.on_changed(slider_changed)

def slider_td_changed(val):
    lines_data = CalcSpectra(td=val, m=slider.val)
    for i,l in enumerate(lines):
      l.set_ydata(lines_data[:,i])
    fig.canvas.draw_idle()
slider_td.on_changed(slider_td_changed)

plot.show()

# Plot the bulk band structure, and compare it to the finite band structure.
import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.widgets import Slider

np.set_printoptions(precision=1,edgeitems=30, linewidth=100000)

a_x = 1   # Lattice constant
#m = -0.1  # Mass term (onsite energy difference)
kx_ = np.arange(-np.pi, np.pi*1.01, np.pi/100)
N = 50 # sites in the y direction

def CalcSpectra(m):
  plotData = np.zeros([len(kx_),2*N])

  for kx_i, kx in enumerate(kx_):
    H = np.zeros([2*N, 2*N])
    for i in range(int(N)):
      for j in range(int(N)):
        if i == j:
          H[i*2][j*2] = m + np.cos(kx*a_x)
          H[i*2+1][j*2+1] = -m - np.cos(kx*a_x)
          H[i*2+1][j*2] = np.sin(kx*a_x)
          H[i*2][j*2+1] = np.sin(kx*a_x)
        if i == j+1:
          H[i*2][j*2] = 0.5
          H[i*2+1][j*2+1] = -0.5
          H[i*2+1][j*2] = -0.5
          H[i*2][j*2+1] = 0.5
        if i+1 == j:
          H[i*2][j*2] = 0.5
          H[i*2+1][j*2+1] = -0.5
          H[i*2+1][j*2] = 0.5
          H[i*2][j*2+1] = -0.5
    energies, eigenstates = np.linalg.eig(H)
    plotData[kx_i] = np.sort(energies)

  return plotData


fig = plot.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3)
lines = ax.plot(kx_,CalcSpectra(0.0))
plot.title("Fritz Model Energies")
plot.xlabel("$k_x$")
plot.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],["$-\\pi$","$-\\frac{\\pi}{2}$","$0$","$\\frac{\\pi}{2}$","$\\pi$"])
plot.ylabel("Energy")
slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slider = Slider(slider_ax, '$m$', -2.5, 2.5, valinit=0.0)



def slider_changed(val):
    m = val
    lines_data = CalcSpectra(m)
    for i,l in enumerate(lines):
      l.set_ydata(lines_data[:,i])
    fig.canvas.draw_idle()
slider.on_changed(slider_changed)

plot.show()
print("hi")
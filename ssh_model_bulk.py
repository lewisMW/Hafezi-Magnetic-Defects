# Plot the bulk band structure, and compare it to the finite band structure.
from tkinter import W
import numpy as np
import matplotlib.pyplot as plot
from matplotlib.widgets import Slider

w = 1 #external   #metal
v = 1             #metal
#w=0.5            #insulator
#v=1              #insulator

k = np.linspace(0,2*np.pi,200)
E_plus = np.sqrt(2*v*w*np.cos(k)+v**2+w**2)
E_minus = -E_plus

fig = plot.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.25)
    
[line1] = ax.plot(k, E_plus)
[line2] = ax.plot(k, E_minus)
plot.ylim([-3,3])
plot.xlabel("k")
plot.ylabel("E")
plot.title("SSH bulk Band Structure, w=%d"%w)
plot.xticks([0,np.pi/2,np.pi,3*np.pi/2,2*np.pi],["0","$\\frac{\\pi}{2}$","$\\pi$","$\\frac{3\\pi}{2}$","$2\\pi$"])

slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slider = Slider(slider_ax, 'v', 0.0, 2.0, valinit=v)

def slider_changed(val):
    v = val
    E_plus = np.sqrt(2*v*w*np.cos(k)+v**2+w**2)
    E_minus = -E_plus
    line1.set_ydata(E_plus)
    line2.set_ydata(E_minus)
    fig.canvas.draw_idle()
slider.on_changed(slider_changed)

plot.show()

print("Gap = 2|w-v| = %f" % (2*abs(w-v)))
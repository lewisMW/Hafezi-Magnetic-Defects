#Plot the eigenstates (wave functions) of the finite SSH model.
import numpy as np
import matplotlib.pyplot as plot

w = 1 # External
N = 80
v_pts = 20
plotM = np.zeros([v_pts,2*N])

i=0
for v in np.linspace(0,2,v_pts):  # Internal

    E = np.array([[0,v],[v,0]])
    T = np.array([[0,0],[w,0]])
    T_ = np.array([[0,w],[0,0]])
    o = np.array([[0,0],[0,0]])

    A = np.block([ [o]*(n-1)+[T_]*(1 if n>0 else 0)+[E]+[T]*(1 if n <N-1 else 0)+[o]*(N-n-2) for n in range(N)])

    energies, eigenstates = np.linalg.eig(A)
    
    indexed_energies = sorted(enumerate(energies), key=lambda x: x[1])
    plot.clf()
    for E in indexed_energies:
        print(E[1])
        plot.plot(eigenstates[:,E[0]])
        plot.show()

    plotM[i] = np.sort(energies)
    i += 1
    
plot.plot(plotM)
plot.title("SSH Model Energies, Sweeping Intercell Parameter v. N=%d, w=%d"%(N,w))
plot.xlabel("Intercell Parameter v")
plot.ylabel("Energy")
plot.show()
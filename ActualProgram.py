import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz

def f(u: np.ndarray, alpha):
    res = np.zeros(u.shape)
    res[0] = u[1]-2*u[0]+alpha*(u[1]-u[0])**2-alpha*u[0]**2
    for i in range(1, u.shape[0]-1):
        res[i] = u[i+1]-2*u[i]+u[i-1]+alpha*(u[i+1]-u[i])**2-alpha*(u[i]-u[i-1])**2

    res[-1] = -2*u[-1]+u[-2]+alpha*u[-1]**2-alpha*(u[-1]-u[-2])**2

    return res




NA = 32

spiny = [0]*(NA-1)
spiny[0:2] = [2, -1]

A = toeplitz(spiny)

eigvals, eigvecs = np.linalg.eig(A)

# sort eigenvalues and vectors
args = np.argsort(eigvals)
eigvals = eigvals[args]
eigvecs = eigvecs[:,args]

# initial conditions
Nt = 50000
delta = np.sqrt(1/8)
T = np.linspace(0, Nt*delta, Nt)
u = 4*np.sqrt(2/NA)*np.sin(np.pi*np.arange(1, NA)/NA)
v = np.zeros(NA-1)
E_0 = np.zeros(Nt)

for i in range(Nt):
    F1 = f(u, 0.25)

    #calculate energy at timestep
    xi_0 = np.dot(eigvecs[0], u)
    xip_0 = np.dot(eigvecs[0], v)
    E_0[i] = 0.5*(xip_0 ** 2 + (xi_0 ** 2) * (eigvals[0] ** 2))

    u = u + v*delta + 0.5*F1*delta**2
    F2 = f(u, 0.25)
    v = v + 0.5*delta*(F2+F1)


plt.plot((T*np.sqrt(eigvals[0]))/(2*np.pi), 100*E_0, label='E0')

plt.xlim(0, 160)

plt.show()





diag_A = np.linalg.inv(eigvecs) @ A @ eigvecs



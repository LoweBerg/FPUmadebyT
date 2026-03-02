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




NA = 5

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
u = np.sqrt(2/NA)*np.sin(np.pi*np.arange(1, NA)/5)
v = np.zeros(NA)

for i in range(Nt):
    F1 = f(u, delta)
    u = u + v*delta + 0.5*F1*delta**2
    F2 = f(u, delta)
    v = v + 0.5*delta*(F2+F1)

print(eigvecs)




diag_A = np.linalg.inv(eigvecs) @ A @ eigvecs



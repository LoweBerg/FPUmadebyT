import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz

def f(u: np.ndarray, alpha):
    res = np.zeros(u.shape, dtype=np.float64)
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

# Transpose vectors for shorter indexing later
eigvecs = eigvecs.T

# sort eigenvalues and vectors
args = np.argsort(eigvals)
eigvals = eigvals[args]
eigvecs = eigvecs[args]

# initial conditions
Nt = 50000
delta = np.sqrt(1/8)
T = np.zeros(Nt, dtype=np.float64)
u = 4*eigvecs[0]
v = np.zeros(NA-1, dtype=np.float64)
E = np.zeros((4, Nt), dtype=np.float64)
a = 0.25

plt.figure()
for i in range(4):
    plt.plot(np.arange(1, NA), eigvecs[i] + 0.5*i)
plt.title("Initial conditions")

for i in range(Nt):
    F1 = f(u, a)

    # wanted to do this with no loops but ran into precision errors (I think)
    for j in range(4):
        # calculate energy at timestep
        xi = np.dot(eigvecs[j], u)
        xip = np.dot(eigvecs[j], v)

        # save energy and timestep
        Energy = 100*0.5*(np.pow(xip, 2) + np.dot(np.pow(xi, 2), eigvals[j]))
        E[j, i] = Energy

    T[i] = i*delta*np.sqrt(eigvals[0])/(2*np.pi)

    u = u + v*delta + 0.5*F1*(delta**2)
    F2 = f(u, a)
    v = v + 0.5*delta*(F2+F1)


plt.figure()
for i in range(4):
    plt.plot(T, E[i], label=f'E{i}')

plt.title("F me in the PUT")

plt.ylim(0, 8)
plt.xlim(0, 160)

plt.legend()

plt.show()

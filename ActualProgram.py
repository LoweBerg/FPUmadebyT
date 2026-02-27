import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz


N = 5

spiny = [0]*(N-1)
spiny[0:2] = [2, -1]

A = toeplitz(spiny)

eigvals, eigvecs = np.linalg.eig(A)

diag_A = np.linalg.inv(eigvecs) @ A @ eigvecs



import numpy as np


# A ska vara 31x31, börja med 5x5 så att allt funk

def makematrix(n, a, b, c):
    matrix = np.zeros((n, n), int)
    np.fill_diagonal(matrix, a)
    np.fill_diagonal(matrix[:-1, 1:], b)
    np.fill_diagonal(matrix[1:, :-1], c)
    return matrix

def F1(u, alpha):
    return u[1] - 2*u[0] + alpha*(u[1]- u[0])**2 - alpha*u[0]**2   # tänk på att index är -1 från formeln här pga python

def Fi(u, alpha, N):
    for i in range(2-1, N-1):     # börjar på i=1 här återigen pga pythons indexering
        u[i+1] - 2*u[i] + u[i-1] + alpha*(u[i+1]-u[i])**2 - alpha*(u[i]-u[i-1])**2
    return ?????????????????????????

def FNminus1(u, alpha, N):
    return -2*u[N-2] + u[N-3] + alpha*u[N-2]**2 - alpha*(u[N-2]-u[N-3])**2

# make the matrix
A = makematrix(4, 2, -1, -1)        

# extract its eigenthingies
eigvals, eigvecs = np.linalg.eig(A)  
# note to self: när man checkar ifall eigvals är rätt, så ska N i sin vara +1 i förhållande till As dimensioner!! A har ju dim N-1

diagonalmat = np.diag(eigvals)

# order eigvals and vecs from low to high
eigvals = np.sort(eigvals)
eigvecs = np.sort(eigvecs)      # blir sorterad så att man får högst summa typ

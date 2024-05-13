"""

Functions to solve linear systems of equations.

Kuo-Chuan Pan
2024.05.05

"""
import numpy as np

def solveLowerTriangular(L,b):
    """
    Solve a linear system with a lower triangular matrix L.

    Arguments:
    L -- a lower triangular matrix
    b -- a vector

    Returns:
    x -- the solution to the linear system
    """
    n  = len(b)
    x  = np.zeros(n)
    bs = np.copy(b) # backup of b

    #  can be numba kernel
    for j in np.arange(n):
        if L[j,j] == 0:
            raise ValueError("Matrix is singular.")
        x[j] = bs[j]/L[j,j]

        for i in np.arange(j+1,n):
            bs[i] -= L[i,j]*x[j]
    
    return x


def solveUpperTriangular(U,b):
    """
    Solve a linear system with an upper triangular matrix U.

    Arguments:
    U -- an upper triangular matrix
    b -- a vector

    Returns:
    x -- the solution to the linear system

    """
    n  = len(b)
    x  = np.zeros(n)
    bs = np.copy(b) # backup of b

    #  can be numba kernel
    for j in np.arange(n-1,-1,-1):
        if U[j,j] == 0:
            raise ValueError("Matrix is singular.")
        x[j] = bs[j]/U[j,j]

        for i in np.arange(j):
            bs[i] -= U[i,j]*x[j]
    
    return x


def lu(A):
    """
    Perform LU decomposition on a square matrix A.

    Arguments:
    A -- a square matrix

    Returns:
    L -- a lower triangular matrix
    U -- an upper triangular matrix

    """
    n  = len(A)
    L  = np.identity(n)
    U  = np.zeros((n,n))
    M  = np.zeros((n,n))
    As = np.copy(A) # backup of A 

    for k in np.arange(n):
        if As[k,k]==0:
            raise ValueError("Matrix is singular.")
        for i in np.arange(k+1,n):
            M[i,k] = As[i,k]/As[k,k]
    
        for j in np.arange(k+1,n):
            for i in np.arange(k+1,n):
                As[i,j] -= M[i,k]*As[k,j]

    # compute L and U
    for i in np.arange(n):
        L[i,:i] = M[i,:i]
        U[i,i:] = As[i,i:]

    return L, U


def lu_solve(A,b):
    """
    Solve a linear system with a square matrix A using LU decomposition.

    Arguments:
    A -- a square matrix
    b -- a vector

    Returns:
    x -- the solution to the linear system

    """

    x = np.zeros(len(b))
    l, u = lu(A)
    # L(U x) = b; U x = y
    # L y = b
    y = solveLowerTriangular(l,b)

    # y = U x
    x = solveUpperTriangular(u,y)
    return x
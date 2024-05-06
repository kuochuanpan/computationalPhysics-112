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
    bs = np.copy(b)

    for j in range(n):
        # check if L[i,i] is singular
        if L[j,j] == 0:
            raise ValueError("L[{},{}] is zero".format(i,i))
        x[j] = bs[j]/L[j,j]
        
        for i in range(j,n):
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
    bs = np.copy(b)

    for j in range(n-1,-1,-1):
        # check if U[i,i] is singular
        if U[j,j] == 0:
            raise ValueError("U[{},{}] is zero".format(i,i))
        x[j] = bs[j]/U[j,j]
        
        for i in range(j):
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
    As = np.copy(A)

    for k in np.arange(n):
        # check if A[i,i] is singular
        if As[k,k] == 0:
            raise ValueError("A[{},{}] is zero".format(i,i))
        
        for i in np.arange(k+1,n):
            M[i,k] = As[i,k]/As[k,k]

        for j in np.arange(k+1,n):
            for i in range(k+1,n):
                As[i,j] -= M[i,k]*As[k,j]

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

    # First, we decompose A into L and U
    # using the LU decomposition

    L, U = lu(A)

    # solve L y = b
    y    = solveLowerTriangular(L,b)

    # solve U x = y
    x    = solveUpperTriangular(U,y)

    return x
#!/usr/bin/python

from __future__ import division, print_function
from math import ceil,log
import numpy as np

# table 10.4
b_d = {
    3:  [120,60,12,1],
    5:  [30240,15120,3360,420,30,1],
    7:  [17297280, 8648640,1995840, 277200, 25200,1512, 56,1],
    9:  [17643225600,8821612800,2075673600,302702400,30270240, 2162160,110880,3960,90,1],
    13: [64764752532480000,32382376266240000,7771770303897600,1187353796428800,129060195264000,10559470521600,670442572800,33522128640,1323241920,40840800,960960,16380,182,1]
}

# table 10.2
theta3  = 1.5e-2
theta5  = 2.5e-1
theta7  = 9.5e-1
theta9  = 2.1e0
theta13 = 5.4e0


def _expm_pade(A,M):
    dtype = A.dtype
    dim,dim = A.shape
    b = b_d[M]

    U = b[1]*np.eye(dim, dtype=dtype)
    V = b[0]*np.eye(dim, dtype=dtype)

    A2 = np.dot(A,A)
    A2n = np.eye(dim, dtype=dtype)

    # evaluate (10.33)
    for i in range(1,M//2+1):
        A2n = np.dot(A2n,A2)
        U += b[2*i+1]*A2n
        V += b[2*i]  *A2n

    del A2,A2n
    U = np.dot(A,U)

    return np.linalg.solve(V-U, V+U)


def _expm_ss(A,norm):
    # algorithm 10.20, from line 7
    dim,dim = A.shape
    b = b_d[13]

    s = max(0, int(ceil(log(norm/theta13)/log(2))))
    if s > 0:
        A /= 2**s

    Id = np.eye(dim)
    A2 = np.dot(A,A)
    A4 = np.dot(A2,A2)
    A6 = np.dot(A2,A4)

    U = np.dot(A, np.dot(A6, b[13]*A6+b[11]*A4+b[9]*A2)+b[7]*A6+b[5]*A4+b[3]*A2+b[1]*Id)
    V = np.dot(A6, b[12]*A6+b[10]*A4+b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*Id

    del A,Id,A2,A4,A6

    r13 = np.linalg.solve(V-U, V+U)

    del V,U

    return np.linalg.matrix_power(r13, 2**s)


def expm(A):
    """
    Calculate the matrix exponential of a square matrix A: MatrixExp(A)

    This module implements algorithm 10.20 from [1]. The matrix exponential is
    calculated using scaling and squaring, and a Pade approximation.

    [1] Functions of Matrices: Theory and Computation, Nicholas J. Higham, 2008
    """
    rows, columns = A.shape
    if rows != columns:
        raise BaseException("A must be a square matrix")

    # calculate the norm of A
    norm = np.linalg.norm(A, ord=1)

    if   norm < theta3:
        return _expm_pade(A,3)
    elif norm < theta5:
        return _expm_pade(A,5)
    elif norm < theta7:
        return _expm_pade(A,7)
    elif norm < theta9:
        return _expm_pade(A,9)
    else:
        return _expm_ss(A,norm)


if __name__ == "__main__":
    from time import time
    from scipy import linalg
    A = np.random.random((2000,2000))/10

    print("norm ", np.linalg.norm(A, ord=1))

    t = time()
    X_scipy = linalg.expm(A)
    t = time()-t
    print("scipy", t)

    print()

    t = time()
    X_expm  = expm(A)
    t = time()-t
    print("expm ", t)
    print()
    print("diff ", np.linalg.norm(X_expm-X_scipy, ord=1))
    print("diff ", np.linalg.norm(X_expm-X_scipy, ord=1)/np.linalg.norm(X_scipy))

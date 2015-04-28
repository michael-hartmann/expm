#!/usr/bin/python

from __future__ import division
import numpy as np
from scipy import linalg

b_d = {
    3:  [120,60,12,1],
    5:  [30240,15120,3360,420,30,1],
    7:  [17297280, 8648640,1995840, 277200, 25200,1512, 56,1],
    9:  [17643225600,8821612800,2075673600,302702400,30270240, 2162160,110880,3960,90,1],
    13: [64764752532480000,32382376266240000,7771770303897600,1187353796428800,129060195264000,10559470521600,670442572800,33522128640,1323241920,40840800,960960,16380,182,1]
}

theta3 = 1.5e-2
theta5 = 2.5e-1
theta7 = 9.5e-1


# M is odd
def _expm_pade(A,M):
    dim,dim = A.shape
    b = b_d[M]

    U = b[1]*np.eye(dim)
    V = b[0]*np.eye(dim)

    # evaluate (10.33)
    for m in range(3,M+1,2):
        print m
        U += b[m]*np.linalg.matrix_power(A,m-1)
    U = np.dot(A,U)

    for m in range(2,M+1,2):
        V += b[m]*np.linalg.matrix_power(A,m)

    VpU = V+U
    VmU = V-U

    inv = np.linalg.inv(VmU)

    return np.dot(inv, VpU)


def expm(A):
    rows, columns = A.shape
    if rows != columns:
        raise BaseException("A must be square matrix")

    # calculate the norm of A
    norm = np.linalg.norm(A)

    if   norm < theta3:
        return _expm_pade(A,3)
    elif norm < theta5:
        return _expm_pade(A,5)
    elif norm < theta7:
        return _expm_pade(A,7)
    else:
        raise BaseException("Not implemented yet")


if __name__ == "__main__":
    A = np.array([[1,2],[3,4]])/6
    print "norm A", np.linalg.norm(A)
    print
    X_scipy = linalg.expm(A)
    X_expm  = expm(A)
    print "scipy", X_scipy
    print
    print "expm ", X_expm
    print
    print "diff ", X_expm-X_scipy, np.linalg.norm(X_expm-X_scipy)

    

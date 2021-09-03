# -*- coding: utf-8 -*-
"""
Linear algebra.
"""

import numpy as np

def get_spectral_radius_and_eigenvector(M):
    """
    Does what it says

    Parameters
    ==========
    M : numpy.ndarray
        Array of square shape

    Returns
    =======
    rho : numpy.float64
        spectral radius
    y : numpy.ndarray
        corresponding eigenvector, normalized such that

        .. math::

            \sum_i y_i = 1
    """

    shape = M.shape
    assert(shape[0]==shape[1])

    rho, y = np.linalg.eig(M)
    rho = np.abs(rho)

    ndx = np.argmax(rho)

    rho = rho[ndx]
    y = np.abs(y[:,ndx])
    y /= sum(y)

    return rho, y

if __name__=="__main__":
    M = np.array([[ 1.0,0.],[0.,2.0]])
    print(get_spectral_radius_and_eigenvector(M))

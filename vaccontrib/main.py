# -*- coding: utf-8 -*-
"""
Load covid data.
"""

import numpy as np
from vaccontrib.linalg import get_spectral_radius_and_eigenvector

def get_next_generation_matrix_from_matrices(R0,gamma, S, N, s, r, a, b):
    """
    """

    M, _ = gamma.shape
    _, V = S.shape

    a0 = a[:,0]
    b0 = b[:,0]

    K0 = gamma.dot(np.diag(a0)).dot(np.diag(1/b0))
    rho0, _ = get_spectral_radius_and_eigenvector(K0)

    if not hasattr(R0,'__len__'):
        R0 = np.ones(V) * R0

    K = 1/rho0 * \
            R0[None,None,None,:] * \
            gamma[:,:,None,None] * \
            (1-s[:,None,:,None]) * \
            (1-r[None,:,None,:]) * \
            S[:,None,:,None] * \
            1 / N[None,:,None,None] * \
            a[None,:,None,:] * \
            1 / b[None,:,None,:]

    return K

def get_4d_matrix_as_2d_block(K):
    M, _, V, __ = K.shape
    _K = np.zeros((M*V, M*V))
    for i in range(M):
        for j in range(M):
            _K[i*V:(i+1)*V,j*V:(j+1)*V] = K[i,j,:,:]
    return _K

def get_contribution_matrix(K):

    M, _, V, __ = K.shape
    _K = get_4d_matrix_as_2d_block(K)
    R, y = get_spectral_radius_and_eigenvector(_K)
    y = y.reshape(M,V)

    C = K * y[None,:,None,:]

    return C

def get_reduced_contribution_matrix(K):
    C = get_contribution_matrix(K)
    C = C.sum(axis=0).sum(axis=0)
    return C

def get_reduced_vaccinated_susceptile_contribution_matrix(K):
    C = get_reduced_contribution_matrix(K)
    _C = np.zeros((2,2))
    _C[0,0] = C[0,0]
    _C[0,1] = C[0,1:].sum()
    _C[1,0] = C[1:,0].sum()
    _C[1,1] = C[1:,1:].sum()
    return _C

def get_reduced_population_contribution_matrix(K):
    C = get_contribution_matrix(K)
    C = C.sum(axis=-1).sum(axis=-1)
    return C

if __name__=="__main__":

    gamma = np.array([[1.,1.],[1.,1.]])
    S = np.array([[0.4,0.6],[0.4,0.6]])
    N = np.array([1.,1.])
    s = np.array([[0.,1],[0.,1]])
    r = np.array([[0.,0.],[0.,0.]])
    a = np.array([[1.,1.],[1.,1.]])
    b = np.array([[1.,1.],[1.,1.]])

    K = get_next_generation_matrix_from_matrices(4,gamma,S,N,s,r,a,b)
    C = get_contribution_matrix(K)
    #print(get_spectral_radius_and_eigenvector(K.reshape(2*2,2*2)))
    print(C)
    #print(C.reshape((4,4)))
    print(C.sum(axis=0).sum(axis=0))
    print(C.sum())


    #print(C[:,:,0,1])

    #print()
    #print(K[:,:,0,1])
    #print(C)

    rows = []
    for row in np.array_split(C,2,axis=0):
        rows.append([])
        for col in np.array_split(row,2,axis=1):
            rows[-1].append(col.reshape(2,2))

    #print(np.array_split(C,2,axis=0))
    #print(np.block(rows))




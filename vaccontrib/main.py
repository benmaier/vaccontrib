# -*- coding: utf-8 -*-
"""
Load covid data.
"""

import numpy as np
from vaccontrib.linalg import get_spectral_radius_and_eigenvector

def get_next_generation_matrix_from_matrices(R0,gamma, S, N, s, r, a, b):
    """
    """

    assert(N.ndim==1)
    for matrix in [gamma,S,s,r,a,b]:
        assert(matrix.ndim==2)

    M, _ = gamma.shape
    _, V = S.shape

    a0 = a[:,0]
    b0 = b[:,0]

    K0 = gamma.dot(np.diag(a0)).dot(np.diag(1/b0))
    rho0, _ = get_spectral_radius_and_eigenvector(K0)

    if not hasattr(R0,'__len__'):
        R0 = np.ones(V) * R0
    else:
        R0 = np.array(R0,dtype=np.float64)

            #np.sqrt(R0[None,None,None,:]*R0[None,None,:,None]) * \
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

def get_homogeneous_contribution_matrix(R0,v,r,s):

    if not hasattr(R0,'__len__'):
        _R = np.ones(2) * R0
    else:
        _R = np.array(R0,dtype=np.float64)

    assert(_R.ndim==1)

    _S = 0
    _V = 1

    R = np.zeros((2,2))

    R[_S,_S] = (1-v)**2/(1-v*s) * _R[_S]
    R[_S,_V] = v*(1-v)*(1-s)*(1-r)/(1-v*s) * _R[_V]
    R[_V,_S] = v*(1-v)*(1-s)/(1-v*s) * _R[_S]
    R[_V,_V] = v**2 * (1-s)**2 * (1-r)/(1-v*s) * _R[_V]

    return R

def get_4d_matrix_as_2d_block(K):
    M, _, V, __ = K.shape
    _K = np.zeros((M*V, M*V))
    for i in range(M):
        for j in range(M):
            _K[i*V:(i+1)*V,j*V:(j+1)*V] = K[i,j,:,:]
    return _K

def get_contribution_matrix(K,return_eigenvector_too=False):

    M, _, V, __ = K.shape
    _K = get_4d_matrix_as_2d_block(K)
    R, y = get_spectral_radius_and_eigenvector(_K)
    y = y.reshape(M,V)

    C = K * y[None,:,None,:]

    if return_eigenvector_too:
        return C, y
    else:
        return C

def get_reduced_contribution_matrix(K):
    C = get_contribution_matrix(K)
    C = C.sum(axis=0).sum(axis=0)
    return C

def get_reduced_vaccinated_susceptible_contribution_matrix(K):
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

def get_eigenvector(K):
    _, y = get_contribution_matrix(K,return_eigenvector_too=True)
    return y

def get_reduced_vaccinated_susceptible_eigenvector(K):
    _, y = get_contribution_matrix(K,return_eigenvector_too=True)
    y = y.sum(axis=0)
    return y

def get_reduced_population_eigenvector(K):
    _, y = get_contribution_matrix(K,return_eigenvector_too=True)
    y = y.sum(axis=-1)
    return y

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

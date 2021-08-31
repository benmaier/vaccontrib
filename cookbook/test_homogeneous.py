from vaccontrib.main import (
        get_next_generation_matrix_from_matrices,
        get_reduced_vaccinated_susceptible_eigenvector,
        get_homogeneous_contribution_matrix,
    )

import numpy as np
_v = 0.69
_s = 0.8
_r = 0.5
R0 = [1.,4.]
R1 = [4.,4.]

gamma = np.array([[1.]])
S = np.array([[1-_v,_v]])
N = np.array([1.])
s = np.array([[0.,_s]])
r = np.array([[0.,_r]])
a = np.array([[1.,1.]])
b = np.array([[1.,1.]])

K0 = get_next_generation_matrix_from_matrices(R0,gamma,S,N,s,r,a,b)
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)

print(y0)

K1 = get_next_generation_matrix_from_matrices(R1,gamma,S,N,s,r,a,b)
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)

print(y1)

from vaccontrib.covid import get_reduced_vaccinated_susceptible_contribution_matrix_covid, get_reduced_population_contribution_matrix_covid

from vaccontrib.covid import get_next_generation_matrix_covid
from vaccontrib.main import get_reduced_vaccinated_susceptible_eigenvector, get_reduced_vaccinated_susceptible_contribution_matrix, get_next_generation_matrix_from_matrices, convert_4d_matrix_to_2d_block


import numpy as np
np.set_printoptions(linewidth=100)
_a = np.array

gamma = _a([[ 1., 1.],
            [ 10., 1.]])
S = _a([ [ 0.3, 0.7] ,
         [ 1.5, 0.5] ])
N = _a([ 1., 2.])
s = _a([ [0., 0.8],
         [0., 0.2] ])
r = _a([ [0., 0.8],
         [0., 0.2] ])
a = _a([ [1., 1],
         [1., 1] ])
b = _a([ [1., 1],
         [1., 1] ])

R0 = 6

K0 = get_next_generation_matrix_from_matrices(R0, gamma, S, N, s, r, a, b)
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)
C0 = get_reduced_vaccinated_susceptible_contribution_matrix(K0)


s[:,0] = 1-np.sqrt(1/6)
r[:,0] = 1-np.sqrt(1/6)
s[:,0] = 1
r[:,0] = 1

K1 = get_next_generation_matrix_from_matrices(R0, gamma, S, N, s, r, a, b)
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)
C1 = get_reduced_vaccinated_susceptible_contribution_matrix(K1)

#print("with R =", R0)
print("K_reduced =\n", K0.sum(axis=0).sum(axis=0))
print("K =\n", convert_4d_matrix_to_2d_block(K0))
print("C =\n", C0)
print("y =\n", y0)
print("R_eff =", C0.sum(axis=0))
print()
#print("with R =", R1)
print("K_reduced =\n", K1.sum(axis=0).sum(axis=0))
print("K =\n", convert_4d_matrix_to_2d_block(K1))
print("C =\n", C1)
print("y =\n", y1)
print("R_eff =", C1.sum(axis=0))
print()

_K = K1.sum(axis=0).sum(axis=0)
_K = convert_4d_matrix_to_2d_block(K1)
print(_K.dot(np.ones(2*2)))

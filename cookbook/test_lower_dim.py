from vaccontrib.covid import get_reduced_vaccinated_susceptible_contribution_matrix_covid, get_reduced_population_contribution_matrix_covid

from vaccontrib.covid import get_next_generation_matrix_covid
from vaccontrib.main import get_reduced_vaccinated_susceptible_eigenvector, get_reduced_vaccinated_susceptible_contribution_matrix, get_next_generation_matrix_from_matrices, get_4d_matrix_as_2d_block


import numpy as np
_a = np.array

gamma = _a([[ 1., 1.],
            [ 2., 1.]])
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

R0 = [6,6]
K0 = get_next_generation_matrix_from_matrices(R0, gamma, S, N, s, r, a, b)
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)
C0 = get_reduced_vaccinated_susceptible_contribution_matrix(K0)


R1 = [1,6]
K1 = get_next_generation_matrix_from_matrices(R1, gamma, S, N, s, r, a, b)
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)
C1 = get_reduced_vaccinated_susceptible_contribution_matrix(K1)

print("with R =", R0)
print("C =\n", C0)
print("y =\n", y0)
print()
print("with R =", R1)
print("C =\n", C1)
print("y =\n", y1)
print()

gamma = _a([[ 1., 1., 2.],
            [ 1., 1., 1.],
            [ 1., 3., 1.],
            ])
S = _a([ [ 0.3, 0.4, 0.3] ,
         [ 1.5, 0.3, 0.2] ,
         [ 2.0, 0.6, 0.4] ,
       ])
N = _a([ 1., 2., 3.])
s = _a([ [0., 0.8, 0.7],
         [0., 0.2, 0.3],
         [0., 0.5, 0.2],
         ])
r = _a([ [0., 0.6, 0.3],
         [0., 0.4, 0.8],
         [0., 0.3, 0.9],
         ])
a = _a([ [0.5, 0.5, 0.5],
         [1., 1, 1],
         [1., 1, 1],
         ])
b = _a([ [1., 1.5, 1.5],
         [1., 1.5, 1.5],
         [1., 1.5, 1.5],
         ])

R0 = [6,6,6]
K0 = get_next_generation_matrix_from_matrices(R0, gamma, S, N, s, r, a, b)
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)
C0 = get_reduced_vaccinated_susceptible_contribution_matrix(K0)


R1 = [1,6,6]
K1 = get_next_generation_matrix_from_matrices(R1, gamma, S, N, s, r, a, b)
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)
C1 = get_reduced_vaccinated_susceptible_contribution_matrix(K1)

print("with R =", R0)
print("K =\n", get_4d_matrix_as_2d_block(K0))
print("C =\n", C0)
print("y =\n", y0)
print()
print("with R =", R1)
print("K =\n", get_4d_matrix_as_2d_block(K1))
print("C =\n", C1)
print("y =\n", y1)
print()

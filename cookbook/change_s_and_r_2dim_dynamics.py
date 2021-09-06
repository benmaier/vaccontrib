from vaccontrib.covid import get_reduced_vaccinated_susceptible_contribution_matrix_covid, get_reduced_population_contribution_matrix_covid

from vaccontrib.covid import get_next_generation_matrix_covid
from vaccontrib.main import get_reduced_vaccinated_susceptible_eigenvector, get_reduced_vaccinated_susceptible_contribution_matrix, get_next_generation_matrix_from_matrices, convert_4d_matrix_to_2d_block

from vaccontrib.linalg import get_spectral_radius_and_eigenvector

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

R0 = 2

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


_K0 = convert_4d_matrix_to_2d_block(K0)
_K1 = convert_4d_matrix_to_2d_block(K1)
_, y0_full = get_spectral_radius_and_eigenvector(_K0)
_, y1_full = get_spectral_radius_and_eigenvector(_K1)

print("K_reduced =\n", K0.sum(axis=0).sum(axis=0))
print("K =\n", convert_4d_matrix_to_2d_block(K0))
print("C =\n", C0)
print("y_reduced =\n", y0)
print("y0_full =\n", y0_full)
print("R_eff =", C0.sum(axis=0))
print("C_normed =\n", C0/C0.sum())
print()
#print("with R =", R1)
print("K_reduced =\n", K1.sum(axis=0).sum(axis=0))
print("K =\n", convert_4d_matrix_to_2d_block(K1))
print("C =\n", C1)
print("y_reduced =\n", y1)
print("y0_full =\n", y1_full)
print("R_eff =", C1.sum(axis=0))
print("C_normed =\n", C1/C1.sum())
print()


y = np.array([1,1,1,1.])
growth = [1.]
VV_growth = [1.]
ytot = sum(y)
ytots = [ytot]
y_d_eig = []

for K in [_K0,_K1]:
    _, eig = get_spectral_radius_and_eigenvector(K)

    for g in range(10):
        y_d_eig.append((y/np.linalg.norm(y)).dot(eig/np.linalg.norm(eig)))
        ytot_old = ytot
        yold = y.copy()
        y = K.dot(y)
        ytot = y.sum()
        ytots.append(ytot)
        growth.append(ytot/ytot_old)
        VV_growth.append(y[[1,3]].sum()/yold[[1,3]].sum())

y_d_eig.append((y/np.linalg.norm(y)).dot(eig/np.linalg.norm(eig)))
import matplotlib.pyplot as pl

fig, ax = pl.subplots(2,2,figsize=(8,4))
ax[0,0].plot(growth)
ax[0,0].plot([0,len(growth)-1],[C1.sum()]*2,':')
ax[0,0].plot([0,len(growth)-1],[C0.sum()]*2,':')

ax[1,0].plot(VV_growth)
ax[1,0].plot([0,len(VV_growth)-1],[C1[1,1]]*2,':')
ax[1,0].plot([0,len(VV_growth)-1],[C0[1,1]]*2,':')

ax[0,1].plot(ytots)
ax[0,0].set_yscale('log')


ax[1,1].plot(y_d_eig)

pl.show()

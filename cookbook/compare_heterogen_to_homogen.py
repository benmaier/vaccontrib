from vaccontrib.covid import get_reduced_vaccinated_susceptible_contribution_matrix_covid, get_reduced_population_contribution_matrix_covid

from vaccontrib.covid import get_next_generation_matrix_covid
from vaccontrib.main import get_reduced_vaccinated_susceptible_eigenvector, get_reduced_vaccinated_susceptible_contribution_matrix

print("====== heterogeneous system with heterogeneous vaccination values ======")

R0 = [6,6,6,6,6]
K0 = get_next_generation_matrix_covid(R0,'delta')
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)
C0 = get_reduced_vaccinated_susceptible_contribution_matrix(K0)


R1 = [0,6,6,6,6]
K1 = get_next_generation_matrix_covid(R1,'delta')
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)
C1 = get_reduced_vaccinated_susceptible_contribution_matrix(K1)

print("with R =", R0)
print("C_VS =")
print(C0)
print("R_eff =", C0.sum(axis=0))
print()
print("with R =", R1)
print("C_VS =")
print(C1)
print("R_eff =", C1.sum(axis=0))
print()

print("====== compare to homogeneous system with homogeneous vaccination values ======")

from vaccontrib.io import get_homogeneous_vaccination_parameters
from vaccontrib.main import get_homogeneous_contribution_matrix

v, s, r = get_homogeneous_vaccination_parameters('delta')

print("v", v)
print("s", s)
print("r", r)

_C0 = get_homogeneous_contribution_matrix(R0[:2],v,r,s)
_C1 = get_homogeneous_contribution_matrix(R1[:2],v,r,s)

print("with R =", R0[:2])
print("C_VS =")
print(_C0)
print("R_eff =", _C0.sum(axis=0))
print()
print("with R =", R1[:2])
print("C_VS =")
print(_C1)
print("R_eff =", _C1.sum(axis=0))
print()


print("====== compare to heterogeneous system with homogeneous vaccination values ======")

from vaccontrib.covid import get_covid_matrices
from vaccontrib.main import get_next_generation_matrix_from_matrices


matrices = get_covid_matrices(variant='delta')
matrices['s'][1:,1:] = s
matrices['r'][1:,1:] = r
matrices['S'][:,0] = matrices['N']*(1-v)
matrices['S'][:,1:] = matrices['N'][:,None]*v / 4

R0 = [6,6,6,6,6]
K0 = get_next_generation_matrix_from_matrices(R0, **matrices)
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)
C0 = get_reduced_vaccinated_susceptible_contribution_matrix(K0)


R1 = [0,6,6,6,6]
K1 = get_next_generation_matrix_from_matrices(R1, **matrices)
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)
C1 = get_reduced_vaccinated_susceptible_contribution_matrix(K1)

print("with R =", R0)
print("C_VS =")
print(C0)
print("R_eff =", C0.sum(axis=0))
print()
print("with R =", R1)
print("C_VS =")
print(C1)
print("R_eff =", C1.sum(axis=0))
print()


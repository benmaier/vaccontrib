from vaccontrib.covid import get_reduced_vaccinated_susceptible_contribution_matrix_covid, get_reduced_population_contribution_matrix_covid

from vaccontrib.covid import get_next_generation_matrix_covid
from vaccontrib.main import get_reduced_vaccinated_susceptible_eigenvector, get_reduced_vaccinated_susceptible_contribution_matrix

R0 = [6,6,6,6,6]
K0 = get_next_generation_matrix_covid(R0,'delta')
y0 = get_reduced_vaccinated_susceptible_eigenvector(K0)
C0 = get_reduced_vaccinated_susceptible_contribution_matrix(K0)


R1 = [1,6,6,6,6]
K1 = get_next_generation_matrix_covid(R1,'delta')
y1 = get_reduced_vaccinated_susceptible_eigenvector(K1)
C1 = get_reduced_vaccinated_susceptible_contribution_matrix(K1)

print("with R =", R0)
print(C0)
print()
print("with R =", R1)
print(C1)
print()

from vaccontrib.io import get_homogeneous_vaccination_parameters
from vaccontrib.main import get_homogeneous_contribution_matrix

v, s, r = get_homogeneous_vaccination_parameters('delta')

print("v", v)
print("s", s)
print("r", r)

_C0 = get_homogeneous_contribution_matrix(R0[:2],v,r,s)
_C1 = get_homogeneous_contribution_matrix(R1[:2],v,r,s)

print("with R =", R0[:2])
print(_C0)
print()
print("with R =", R1[:2])
print(_C1)
print()


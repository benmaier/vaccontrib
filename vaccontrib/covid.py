# -*- coding: utf-8 -*-
"""
Load covid data.
"""

import numpy as np
import vaccontrib.io as io

from vaccontrib.linalg import get_spectral_radius_and_eigenvector
from vaccontrib import (
            get_next_generation_matrix_from_matrices,
            get_contribution_matrix,
            get_reduced_contribution_matrix,
            get_reduced_vaccinated_susceptile_contribution_matrix,
            get_reduced_population_contribution_matrix,
        )

def get_next_generation_matrix_covid(R0,variant='alpha'):
    """
    """
    gamma = io.get_contact_matrix()
    S = io.get_disease_free_state()
    N = io.get_population_sizes()
    s = io.get_susceptibility_reduction(variant=variant)
    r = io.get_transmissibility_reduction(variant=variant)
    a = io.get_relative_infection_rate(variant=variant)
    b = io.get_relative_recovery_rate(variant=variant)

    K = get_next_generation_matrix_from_matrices(R0, gamma, S, N, s, r, a, b)

    return K

def get_contribution_matrix_covid(R0,variant='alpha'):
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_contribution_matrix(K)
    return C

def get_reduced_contribution_matrix_covid(R0,variant='alpha'):
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_reduced_contribution_matrix(K)
    return C

def get_reduced_vaccinated_susceptile_contribution_matrix_covid(R0,variant='alpha'):
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_reduced_vaccinated_susceptile_contribution_matrix(K)
    return C

def get_reduced_population_contribution_matrix_covid(R0,variant='alpha'):
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_reduced_population_contribution_matrix(K)
    return C

if __name__=="__main__":

    R0 = np.array([4,4,4,4,4.])
    R0 = np.array([4,4,4,4,4.])
    K1 = get_next_generation_matrix_covid(R0,variant='delta')

    M, _, V, __ = K1.shape

    C = get_contribution_matrix(K1)

    print(C.sum())

    print()
    print()
    print()
    print(get_reduced_vaccinated_susceptile_contribution_matrix_covid(R0,variant='delta'))



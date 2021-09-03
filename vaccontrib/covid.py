# -*- coding: utf-8 -*-
"""
Functions handling covid data.
"""

import numpy as np
import vaccontrib.io as io

from vaccontrib.linalg import get_spectral_radius_and_eigenvector
from vaccontrib import (
            get_next_generation_matrix_from_matrices,
            get_contribution_matrix,
            get_reduced_contribution_matrix,
            get_reduced_vaccinated_susceptible_contribution_matrix,
            get_reduced_population_contribution_matrix,
        )

def get_covid_matrices(variant='alpha'):
    """
    Load all relevant matrices regarding COVID-19
    vaccine efficacy from package data.

    Parameters
    ==========
    variant : str, default = "alpha"
        The variant for which to load data
        ('alpha' or 'delta)

    Returns
    =======
    matrices : dict
        contains

        - ``'gamma'`` : contact matrix, Entry ``contact_matrix[i,j]``
          contains the average number
          of contacts an average `j`-individual has towards
          `i`-individuals.
        - ``'S'`` : disease-free state, Entry ``S[m,v]`` contains
          the number of m-group individuals
          that are in vaccination state ``v``.
        - ``'N'`` :  population sizes,
          Entry ``population_size[m]`` contains the
          size of population `m`.
        - ``'s'`` : susceptibility reduction,
          Entry ``susceptibility_reduction[m,v]`` contains the
          relative susceptibility reduction of individuals of
          vaccination status `v` and population group `m`.
        - ``'r'`` : transmissibility reduction,
          Entry ``transmissibility_reduction[m,v]`` contains the
          relative transmissibility reduction of individuals of
        - ``'a'`` :
          Entry ``relative_infection_rate[m,v]`` contains the
          infection rate (think: shedding rate) of individuals of
          vaccination status `v` and population group `m` relative
          to some base rate.
        - ``'b'`` : relative_recovery rate,
          Entry ``relative_recovery_rate[m,v]`` contains the
          recovery rate of individuals of
          vaccination status `v` and population group `m` relative
          to some base rate.
    """

    gamma = io.get_contact_matrix()
    S = io.get_disease_free_state()
    N = io.get_population_sizes()
    s = io.get_susceptibility_reduction(variant=variant)
    r = io.get_transmissibility_reduction(variant=variant)
    a = io.get_relative_infection_rate(variant=variant)
    b = io.get_relative_recovery_rate(variant=variant)

    return {
                'gamma' : gamma,
                'S' : S,
                'N' : N,
                's' : s,
                'r' : r,
                'a' : a,
                'b' : b,
            }


def get_next_generation_matrix_covid(R0,variant='alpha'):
    """
    Get the next generation matrix of a covid variant.

    Parameters
    ==========
    R0 : float or list of float
        either global reproduction number or a vector
        of reproduction number values, each for in-
        dividuals of a different vaccination status.
    variant : string, default = 'alpha'
        load data for a covid variant from the package data

    Returns
    =======
    K : numpy.ndarray of shape ``M, M, V, V``
        Next generation matrix.
        Entry ``K[i,j,v,w]`` contains the average `(i,v)`-offspring
        of a single `(j,w)`-individual.
    """
    matrices = get_covid_matrices(variant)
    K = get_next_generation_matrix_from_matrices(R0, **matrices)

    return K

def get_contribution_matrix_covid(R0,variant='alpha'):
    """
    Get the contribution matrix of a covid variant.

    Parameters
    ==========
    R0 : float or list of float
        either global reproduction number or a vector
        of reproduction number values, each for in-
        dividuals of a different vaccination status.
    variant : string, default = 'alpha'
        load data for a covid variant from the package data

    Returns
    =======
    C : numpy.ndarray of shape ``M, M, V, V``
        The system's contribution matrix.
        Entry ``C[i,j,v,w]`` contains the average number
        of `(j,w)`-induced `(i,v)`-offspring during exponential
        growth / decay.
    """
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_contribution_matrix(K)
    return C

def get_reduced_contribution_matrix_covid(R0,variant='alpha'):
    """
    Get the reduced contribution matrix of a covid variant
    (where populations were summed out and only
    vaccination stati remain).

    Parameters
    ==========
    R0 : float or list of float
        either global reproduction number or a vector
        of reproduction number values, each for in-
        dividuals of a different vaccination status.
    variant : string, default = 'alpha'
        load data for a covid variant from the package data

    Returns
    =======
    C : numpy.ndarray of shape ``V, V``
        The system's contribution matrix.
        Entry ``C[v,w]`` contains the average number
        of `w`-induced `v`-offspring during exponential
        growth / decay.
    """
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_reduced_contribution_matrix(K)
    return C

def get_reduced_vaccinated_susceptible_contribution_matrix_covid(R0,variant='alpha'):
    """
    Get the reduced contribution matrix of a covid variant
    where populations were summed over and active vaccination
    statuses where summed over, as well, such that only vaccinated/
    not vaccinated remains.

    Parameters
    ==========
    R0 : float or list of float
        either global reproduction number or a vector
        of reproduction number values, each for in-
        dividuals of a different vaccination status.
    variant : string, default = 'alpha'
        load data for a covid variant from the package data

    Returns
    =======
    C : numpy.ndarray of shape ``2, 2``
        The system's contribution matrix.
        Entry ``C[v,w]`` contains the average number
        of `w`-induced `v`-offspring during exponential
        growth / decay where `w` and `v` can be either
        'vaccinated' or 'not vaccinated'
    """
    K = get_next_generation_matrix_covid(R0,variant)
    C = get_reduced_vaccinated_susceptible_contribution_matrix(K)
    return C

def get_reduced_population_contribution_matrix_covid(R0,variant='alpha'):
    """
    Get the reduced contribution matrix of a covid variant
    (where vaccination stati were summed out and only
    population groups remain).

    Parameters
    ==========
    R0 : float or list of float
        either global reproduction number or a vector
        of reproduction number values, each for in-
        dividuals of a different vaccination status.
    variant : string, default = 'alpha'
        load data for a covid variant from the package data

    Returns
    =======
    C : numpy.ndarray of shape ``M, M``
        The system's contribution matrix.
        Entry ``C[i,j]`` contains the average number
        of `i`-induced `j`-offspring during exponential
        growth / decay.
    """
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
    print(get_reduced_vaccinated_susceptible_contribution_matrix_covid(R0,variant='delta'))



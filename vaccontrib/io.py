# -*- coding: utf-8 -*-
"""
Load data regarding COVID-19 vaccine efficacies.
"""

import csv
import numpy as np

from vaccontrib.paths import get_data_dir

_POPULATIONS = ('[00;12)','[12;18)','[18;60)','[60;oo)')
_VACC_STATUSES = ('no','astra','biontech','moderna','jj')

def _array_from_dict(rows,populations,vaccination_statuses):
    """Convert a list of dictionaries to a numpy array (matrix)"""

    M = len(populations)
    V = len(vaccination_statuses)
    data = np.zeros((M, V))
    for ipop, pop in enumerate(populations):
        for ivacc, vacc in enumerate(vaccination_statuses):
            data[ipop,ivacc] = float(rows[pop][vacc])

    return data.squeeze()

def _dict_from_csvfile(fn):
    """Read a CSV-file and convert its content to a list of dictionaries"""

    rows = {}
    with open(fn,'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            age = row.pop('ages')
            rows[age] = row

    return rows

def _get_pop_vacc_csv(fn,
                      populations=_POPULATIONS,
                      vaccination_statuses=_VACC_STATUSES,
                      ):
    """Convert CSV-file data to a numpy matrix"""

    rows = _dict_from_csvfile(fn)
    data = _array_from_dict(rows,populations,vaccination_statuses)

    return data


def get_susceptibility_reduction(fn=None,
                                 variant='alpha',
                                 populations=_POPULATIONS,
                                 vaccination_statuses=_VACC_STATUSES,
                                 ):
    """
    Load susceptibility reduction values from a file. File must have
    structure like this.

    .. code::

        ages,no,astra,biontech,moderna,jj
        [00;12),0,0.73,0.92,0.92,0.73
        [12;18),0,0.73,0.92,0.92,0.73
        [18;60),0,0.73,0.92,0.92,0.73
        [60;oo),0,0.73,0.92,0.92,0.73

    Parameters
    ==========
    fn : string, default = None
        filename to load
    variant : string, default = 'alpha'
        loads data for a covid variant from the package data, if
        ``fn is None``.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    vaccination_statuses : list of string
        list of strings describing the vaccination statuses (first row
        of the csv-data without first column)

    Returns
    =======
    susceptibility_reduction : numpy.ndarray of shape ``M x V``
        Entry ``susceptibility_reduction[m,v]`` contains the
        relative susceptibility reduction of individuals of
        vaccination status `v` and population group `m`.
    """



    if fn is None:
        fn = get_data_dir() / f'susceptibility_reduction_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_transmissibility_reduction(fn=None,
                                 variant='alpha',
                                 populations=_POPULATIONS,
                                 vaccination_statuses=_VACC_STATUSES,
                                 ):
    """
    Load transmissibility reduction values from a file. File must have
    structure like this.

    .. code::

        ages,no,astra,biontech,moderna,jj
        [00;12),0,0.73,0.92,0.92,0.73
        [12;18),0,0.73,0.92,0.92,0.73
        [18;60),0,0.73,0.92,0.92,0.73
        [60;oo),0,0.73,0.92,0.92,0.73

    Parameters
    ==========
    fn : string, default = None
        filename to load
    variant : string, default = 'alpha'
        loads data for a covid variant from the package data, if
        ``fn is None``.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    vaccination_statuses : list of string
        list of strings describing the vaccination statuses (first row
        of the csv-data without first column)

    Returns
    =======
    transmissibility_reduction : numpy.ndarray of shape ``M x V``
        Entry ``transmissibility_reduction[m,v]`` contains the
        relative transmissibility reduction of individuals of
        vaccination status `v` and population group `m`.
    """

    if fn is None:
        fn = get_data_dir() / f'transmissibility_reduction_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_relative_recovery_rate(fn=None,
                               variant='alpha',
                               populations=_POPULATIONS,
                               vaccination_statuses=_VACC_STATUSES,
                              ):
    """
    Load relative recovery rate values from a file. File must have
    structure like this.

    .. code::

        ages,no,astra,biontech,moderna,jj
        [00;12),1,0.73,0.92,0.92,0.73
        [12;18),1,0.73,0.92,0.92,0.73
        [18;60),1,0.73,0.92,0.92,0.73
        [60;oo),1,0.73,0.92,0.92,0.73

    Parameters
    ==========
    fn : string, default = None
        filename to load
    variant : string, default = 'alpha'
        loads data for a covid variant from the package data, if
        ``fn is None``.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    vaccination_statuses : list of string
        list of strings describing the vaccination statuses (first row
        of the csv-data without first column)

    Returns
    =======
    relative_recovery_rate : numpy.ndarray of shape ``M x V``
        Entry ``relative_recovery_rate[m,v]`` contains the
        recovery rate of individuals of
        vaccination status `v` and population group `m` relative
        to some base rate.
    """

    if fn is None:
        fn = get_data_dir() / f'relative_recovery_rate_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_relative_infection_rate(fn=None,
                                variant='alpha',
                                populations=_POPULATIONS,
                                vaccination_statuses=_VACC_STATUSES,
                                ):
    """
    Load relative infection rate values from a file. File must have
    structure like this.

    .. code::

        ages,no,astra,biontech,moderna,jj
        [00;12),1,0.73,0.92,0.92,0.73
        [12;18),1,0.73,0.92,0.92,0.73
        [18;60),1,0.73,0.92,0.92,0.73
        [60;oo),1,0.73,0.92,0.92,0.73

    Parameters
    ==========
    fn : string, default = None
        filename to load
    variant : string, default = 'alpha'
        loads data for a covid variant from the package data, if
        ``fn is None``.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    vaccination_statuses : list of string
        list of strings describing the vaccination statuses (first row
        of the csv-data without first column)

    Returns
    =======
    relative_infection_rate : numpy.ndarray of shape ``M x V``
        Entry ``relative_infection_rate[m,v]`` contains the
        infection rate (think: shedding rate) of individuals of
        vaccination status `v` and population group `m` relative
        to some base rate.
    """

    if fn is None:
        fn = get_data_dir() / f'relative_infection_rate_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_population_sizes(fn=None,
                         populations=_POPULATIONS,
                         header=('number',),
                        ):
    """
    Load sizes of single populations groups from a file.
    File must have structure like this.

    .. code::

        ages,number
        [00;12),9137232
        [12;18),5339517
        [18;60),46495023
        [60;oo),20275029

    Parameters
    ==========
    fn : string, default = None
        filename to load.
        fFnction loads data from the package data, if
        ``fn is None``.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    header : tuple of string, default = ('number',)
        tuple of a single string containing the header keyword
        that points to the column containing population sizes

    Returns
    =======
    population_size : numpy.ndarray of shape ``M``
        Entry ``population_size[m]`` contains the
        size of population `m`.
    """

    if fn is None:
        fn = get_data_dir() / 'population.csv'

    return _get_pop_vacc_csv(fn, populations, header)

def get_fraction_vaccinated(fn=None,
                            populations=_POPULATIONS,
                            header=('fraction_vaccinated',),
                           ):
    """
    Load the fraction of vaccinated individuals per
    population group from a file. File must have
    structure like this.

    .. code::

        ages,fraction_vaccinated
        [00;12),0.
        [12;18),0.206
        [18;60),0.642
        [60;oo),0.838

    Parameters
    ==========
    fn : string, default = None
        filename to load.
        fFnction loads data from the package data, if
        ``fn is None``.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    header : tuple of string, default = ('fraction_vaccinated',)
        tuple of a single string containing the header keyword
        that points to the column containing the fractions vaccinated.

    Returns
    =======
    fraction_vaccinated : numpy.ndarray of shape ``M``
        Entry ``fraction_vaccinated[m]`` contains the relative
        size of the group off vaccinated people in this population
        group.
    """

    if fn is None:
        fn = get_data_dir() / 'vaccinated.csv'

    return _get_pop_vacc_csv(fn, populations, header)

def get_vaccine_fractions(fn=None,
                          populations=_POPULATIONS,
                          header=_VACC_STATUSES[1:],
                         ):
    """
    Load the fractions of vaccines statuses per population
    group from a file.
    File must have structure like this.

    .. code::

        ages,biontech,astra,moderna,jj
        [00;12),0,0,0,0
        [12;18),1,0,0,0
        [18;60),0.69125,0.12925,0.10925,0.07025
        [60;oo),0.69125,0.12925,0.10925,0.07025

    Parameters
    ==========
    fn : string, default = None
        filename to load
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    header : list of string
        list of strings describing the vaccination statuses (first row
        of the csv-data without first column) without "not vaccinated"

    Returns
    =======
    vaccine_fractions : numpy.ndarray of shape ``M x (V-1)``
        Entry ``vaccine_fractions[m,v]`` contains the
        relative fraction of administered full doses of of
        vaccination status `v+1` (relative proportion in population group `m`)
    """

    if fn is None:
        fn = get_data_dir() / 'vaccine_fractions.csv'

    return _get_pop_vacc_csv(fn, populations, header)

def get_contact_matrix(fn=None,
                       populations=_POPULATIONS,
                      ):
    """
    Load the contact matrix from a file.
    File must have structure like this.

    .. code::

        ages,[00;12),[12;18),[18;60),[60;oo)
        [00;12),2.8394495,0.5205262,3.235192,0.6269835
        [12;18),0.8907488,4.4044118,4.745159,0.4811966
        [18;60),0.6357820,0.5449370,6.430791,1.0125184
        [60;oo),0.2825591,0.1267252,2.321924,2.1267606

    Parameters
    ==========
    fn : string, default = None
        filename to load
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)

    Returns
    =======
    contact_matrix : numpy.ndarray of shape ``M x V``
        Entry ``contact_matrix[i,j]`` contains the average number
        of contacts an average `j`-individual has towards
        `i`-individuals.
    """

    if fn is None:
        fn = get_data_dir() / 'contact_matrix.csv'

    return _get_pop_vacc_csv(fn, populations, populations)

def get_disease_free_state(
                            fraction_vaccinated=None,
                            vaccine_fractions=None,
                            population=None,
                            populations=_POPULATIONS,
                            vaccination_statuses=_VACC_STATUSES,
                          ):
    """
    Get the disease-free state of a system defined by vaccine
    fractions, the fraction of vaccinated individuals per
    population group, and the population sizes

    Parameters
    ==========
    fraction_vaccinated : numpy.ndarray of shape ``M``
        Entry ``fraction_vaccinated[m]`` contains the relative
        size of the group off vaccinated people in this population
        group.
    vaccine_fractions : numpy.ndarray of shape ``M x (V-1)``
        Entry ``vaccine_fractions[m,v]`` contains the
        relative fraction of administered full doses of of
        vaccination status `v+1` (relative proportion in population group `m`)
    population : numpy.ndarray of shape ``M``
        Entry ``population[m]`` contains the
        size of population `m`.
    populations : list of string
        list of strings describing the populations (first column
        of the csv-data without first row)
    vaccination_statuses : list of string
        list of strings describing the vaccination statuses (first row
        of the csv-data without first column)

    Returns
    =======
    S : numpy.ndarray of shape ``M x V``
        Entry ``S[m,v]`` contains the number of m-group individuals
        that are in vaccination state ``v``.
    """

    if fraction_vaccinated is None:
        fraction_vaccinated = get_fraction_vaccinated()
    if vaccine_fractions is None:
        vaccine_fractions = get_vaccine_fractions()
    if population is None:
        population = get_population_sizes()

    M = len(populations)
    V = len(vaccination_statuses)

    S = np.zeros((M, V))
    S[:,0] = population * (1-fraction_vaccinated)
    S[:,1:] = (population * fraction_vaccinated)[:,None] * vaccine_fractions

    return S

def get_default_populations():
    """Get the default list of strings defining population groups"""
    return _POPULATIONS

def get_default_vaccination_statuses():
    """Get the default list of strings defining vaccination stati"""
    return _VACC_STATUSES

def get_homogeneous_vaccination_parameters(
                                fraction_vaccinated=None,
                                vaccine_fractions=None,
                                population=None,
                                susceptibility_reduction=None,
                                transmissibility_reduction=None,
                                variant='alpha',
                            ):
    """
    Get the disease-free state of a system defined by vaccine
    fractions, the fraction of vaccinated individuals per
    population group, the population sizes, the susceptibility
    reductions, and the transmissibility reductions.

    Parameters
    ==========
    fraction_vaccinated : numpy.ndarray of shape ``M``
        Entry ``fraction_vaccinated[m]`` contains the relative
        size of the group off vaccinated people in this population
        group.
    vaccine_fractions : numpy.ndarray of shape ``M x (V-1)``
        Entry ``vaccine_fractions[m,v]`` contains the
        relative fraction of administered full doses of of
        vaccination status `v+1` (relative proportion in population group `m`)
    population : numpy.ndarray of shape ``M``
        Entry ``population[m]`` contains the
        size of population `m`.
    susceptibility_reduction : numpy.ndarray of shape ``M x V``
        Entry ``susceptibility_reduction[m,v]`` contains the
        relative susceptibility reduction of individuals of
        vaccination status `v` and population group `m`.
    transmissibility_reduction : numpy.ndarray of shape ``M x V``
        Entry ``transmissibility_reduction[m,v]`` contains the
        relative transmissibility reduction of individuals of
        vaccination status `v` and population group `m`.

    Returns
    =======
    v : float
        Average fraction of individuals that are vaccinated
    s : float
        Average susceptibility reduction of all administered
        full doses
    r : float
        Average transmissibility reduction of all administered
        full doses
    """
    if fraction_vaccinated is None:
        fraction_vaccinated = get_fraction_vaccinated()
    if vaccine_fractions is None:
        vaccine_fractions = get_vaccine_fractions()
    if population is None:
        population = get_population_sizes()

    if susceptibility_reduction is None:
        s = get_susceptibility_reduction(variant=variant)[1:,1:]
    if transmissibility_reduction is None:
        r = get_transmissibility_reduction(variant=variant)[1:,1:]

    M, V = s.shape

    S = np.zeros((M, V))
    S[:,0] = population * (1-fraction_vaccinated)
    S[:,1:] = (population * fraction_vaccinated)[:,None] * vaccine_fractions
    S_vacc = S[1:,1:]

    _v = sum(population * fraction_vaccinated / sum(population))
    _s = (s*S_vacc/S_vacc.sum()).sum()
    _r = (r*S_vacc/S_vacc.sum()).sum()

    return _v, _s, _r

if __name__=="__main__":
    functions = [
                    get_contact_matrix,
                    get_vaccine_fractions,
                    get_fraction_vaccinated,
                    get_population_sizes,
                    get_susceptibility_reduction,
                    get_transmissibility_reduction,
                    get_relative_infection_rate,
                    get_relative_recovery_rate,
                    get_disease_free_state,
                ]

    for f in functions:
        print()
        print(f.__name__)
        print(f())

    S = get_disease_free_state()
    print(S)
    print(S.sum())


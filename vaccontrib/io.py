# -*- coding: utf-8 -*-
"""
Load data.
"""

import csv
import numpy as np

from vaccontrib.paths import get_data_dir

_POPULATIONS = ('[00;12)','[12;18)','[18;60)','[60;oo)')
_VACC_STATUSES = ('no','astra','biontech','moderna','jj')

def _array_from_dict(rows,populations,vaccination_statuses):

    M = len(populations)
    V = len(vaccination_statuses)
    data = np.zeros((M, V))
    for ipop, pop in enumerate(populations):
        for ivacc, vacc in enumerate(vaccination_statuses):
            data[ipop,ivacc] = float(rows[pop][vacc])

    return data.squeeze()

def _dict_from_csvfile(fn):

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

    rows = _dict_from_csvfile(fn)
    data = _array_from_dict(rows,populations,vaccination_statuses)

    return data


def get_susceptibility_reduction(fn=None,
                                 variant='alpha',
                                 populations=_POPULATIONS,
                                 vaccination_statuses=_VACC_STATUSES,
                                 ):

    if fn is None:
        fn = get_data_dir() / f'susceptibility_reduction_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_transmissibility_reduction(fn=None,
                                 variant='alpha',
                                 populations=_POPULATIONS,
                                 vaccination_statuses=_VACC_STATUSES,
                                 ):

    if fn is None:
        fn = get_data_dir() / f'transmissibility_reduction_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_relative_recovery_rate(fn=None,
                               variant='alpha',
                               populations=_POPULATIONS,
                               vaccination_statuses=_VACC_STATUSES,
                              ):

    if fn is None:
        fn = get_data_dir() / f'relative_recovery_rate_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_relative_infection_rate(fn=None,
                                variant='alpha',
                                populations=_POPULATIONS,
                                vaccination_statuses=_VACC_STATUSES,
                                ):

    if fn is None:
        fn = get_data_dir() / f'relative_infection_rate_{variant}.csv'

    return _get_pop_vacc_csv(fn, populations, vaccination_statuses)

def get_population_sizes(fn=None,
                         populations=_POPULATIONS,
                         header=('number',),
                        ):

    if fn is None:
        fn = get_data_dir() / 'population.csv'

    return _get_pop_vacc_csv(fn, populations, header)

def get_fraction_vaccinated(fn=None,
                            populations=_POPULATIONS,
                            header=('fraction_vaccinated',),
                           ):

    if fn is None:
        fn = get_data_dir() / 'vaccinated.csv'

    return _get_pop_vacc_csv(fn, populations, header)

def get_vaccine_fractions(fn=None,
                          populations=_POPULATIONS,
                          header=_VACC_STATUSES[1:],
                         ):

    if fn is None:
        fn = get_data_dir() / 'vaccine_fractions.csv'

    return _get_pop_vacc_csv(fn, populations, header)

def get_contact_matrix(fn=None,
                       populations=_POPULATIONS,
                      ):

    if fn is None:
        fn = get_data_dir() / 'contact_matrix.csv'

    return _get_pop_vacc_csv(fn, populations, populations)

def get_disease_free_state(
                            populations=_POPULATIONS,
                            vaccination_statuses=_VACC_STATUSES,
                          ):

    fraction_vaccinated = get_fraction_vaccinated()
    vaccine_fractions = get_vaccine_fractions()
    population = get_population_sizes()

    M = len(populations)
    V = len(vaccination_statuses)

    S = np.zeros((M, V))
    S[:,0] = population * (1-fraction_vaccinated)
    S[:,1:] = (population * fraction_vaccinated)[:,None] * vaccine_fractions

    return S


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


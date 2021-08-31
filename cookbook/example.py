from vaccontrib.io import (
                    get_contact_matrix,
                    get_vaccine_fractions,
                    get_fraction_vaccinated,
                    get_population_sizes,
                    get_susceptibility_reduction,
                    get_transmissibility_reduction,
                    get_relative_infection_rate,
                    get_relative_recovery_rate,
                    get_disease_free_state,
                )

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


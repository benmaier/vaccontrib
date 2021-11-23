import vaccontrib.io as io

from vaccontrib.covid import (
        get_reduced_vaccinated_susceptible_contribution_matrix_covid,
        get_homogeneous_contribution_matrix_covid,
        get_reduced_eigenvector_covid,
        get_reduced_population_eigenvector_covid,
        get_covid_matrices,
        get_eigenvector_covid,
        get_contribution_matrix_covid,
    )

from vaccontrib.illustration import (
        CircleCapSegmentPresentation,
        get_circular_vector_and_matrix_presentation,
    )

from vaccontrib.main import (
            convert_4d_matrix_to_2d_block
        )
from numpyarray_to_latex import to_ltx

import numpy as np
import matplotlib.pyplot as pl

import argparse

R0 = [3,3]

datavacc = np.array([
            [ 23_119, 159_647, 37_735 ],
            [ 1_102,  66_396, 22_973 ],
        ],dtype=float)

datay = datavacc.copy().T
datay[:,0] -= datay[:,1]
datay = datay / datay.sum()

incidence_distribution_ages = [0.14710685, 0.09388042, 0.60214432, 0.15686842]

if __name__ == "__main__":
    print(f"During the observation period, {datavacc[1,:].sum()/datavacc[0,:].sum()*100}% of symptomatic cases were breakthrough infections")
    print(f"The relative frequency of breakthrough infections per age group were {datavacc[1,:]/datavacc[0,:]=} for age groups [12,18), [18, 60), and 60+")
    print(f"The age distribution of new cases were 14.7%, 9.4%, 60.2%, and 15.7% for age groups [0,12), [12,18), [18, 60), and 60+")
    print(f"The normalized symptomatic COVID-19 population eigenvector, ignoring children and considerung vaccinated and unvaccinated populations, is given as")
    print(datay)

    print()

    parser = argparse.ArgumentParser(description='Compute contribution matrices.')
    parser.add_argument('directories', metavar='dirs', type=str, nargs='+',
                        help='directories for which contributions matrices should be computed')
    parser.add_argument('-u', '--Ru', dest='Ru', type=float,
                        help='Base R-value of unvaccinated',default=1)
    parser.add_argument('-v', '--Rv', dest='Rv', type=float,
                        help='Base R-value of vaccinated',default=1)
    parser.add_argument('-f', '--save-figures', dest='save_figures', action='store_true',
                        help='create, show, and save illustrations',default=False)

    args = parser.parse_args()
    R0 = [args.Ru, args.Rv]

VACC = ('no', 'vacc')



def make_analysis(dirs,Ru, Rv, save_figures = False, R0=1.2, verbose=True):


    resulting_Rs = []
    for data_dir in dirs:

        Rs = np.array([Ru,Rv])
        data_dir = data_dir.rstrip('/')
        matrices = get_covid_matrices('delta',data_dir,vaccination_statuses=VACC)
        N = matrices['N']
        C = get_reduced_vaccinated_susceptible_contribution_matrix_covid(1.,'delta',data_dir,vaccination_statuses=VACC)
        Rs = Rs/C.sum() * R0
        C = get_reduced_vaccinated_susceptible_contribution_matrix_covid(Rs,'delta',data_dir,vaccination_statuses=VACC)
        yred = get_reduced_eigenvector_covid(Rs, 'delta', data_dir,vaccination_statuses=VACC)
        yred = [yred[0], yred[1:].sum()]
        R = C.sum()
        if verbose:
            print("======",data_dir,"======")
            print("Reff =",R)
            print("2x2 relative contribution matrix =\n",C/R)
            print("Contributions made by (u, v)=", C.sum(axis=0)/R)
            print("total proportion of breakthrough infections y[u,v]=", yred)

        y = get_eigenvector_covid(Rs,variant='delta',data_dir=data_dir,vaccination_statuses=VACC)

        # summ vacc status
        newy = np.zeros((3,2))
        # ignore young children
        newy[:,0] = y[1:,0]
        newy[:,1] = y[1:,1:].sum(axis=1)

        proportion_breakthrough_eligible = newy.sum(axis=0)/newy.sum()
        relative_frequency_breakthrough_eligible = newy/newy.sum(axis=1)[:,None]

        if verbose:
            print(f"relative frequency of breakthrough infections in eligible age groups =\n{proportion_breakthrough_eligible}")
            print(f"proportion of breakthrough infections of those eligible to receive vaccine =\n{relative_frequency_breakthrough_eligible}")

            print()
            print(f"{newy/newy.sum()=}")
            print(f"{datay=}")
            print(f"{y.sum(axis=1)=}")


        if save_figures:
            pres = get_circular_vector_and_matrix_presentation(yred, C)
            ax = pres.plot(figsize=(6,6))
            ax.axis('off')
            pres.add_arrows_to_plot()
            pres.add_text_to_plot()

            caps = CircleCapSegmentPresentation(C)
            caps.compute()
            ax = caps.plot()
            ax.axis('off')
            ax.get_figure().savefig(f'caps_{data_dir}.pdf')

        C_full = get_contribution_matrix_covid(Rs,'delta',data_dir,vaccination_statuses=VACC)

        if verbose:
            print()
            print("The full contribution matrix:")
            print(to_ltx(convert_4d_matrix_to_2d_block(C_full/C_full.sum()),
                         fmt='{:4.2f}',
                         separate_columns = range(2,26,2),
                         separate_rows = range(2,26,2),
                         )
                  )
        resulting_Rs.append(R)

    return resulting_Rs

if __name__=="__main__":
    make_analysis(args.directories, args.Ru, args.Rv,args.save_figures)

    pl.show()

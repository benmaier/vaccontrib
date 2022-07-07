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

R0 = [1,1,1,1,1,]
R0 = [3,3]

total = np.array([
            [ 0, 23_119, 159_647, 37_735 ],
            [ 0,  1_102,  66_396, 22_973 ],
        ],dtype=float)

datay = total.T
datay = datay[1:]
datay[:,0] -= datay[:,1]
print(datay.sum())
datay = datay / datay.sum()

print(f"{total[1,:].sum()/total[0,:].sum()=}")
print(f"{total[0,1:].sum()=}")
print(f"{total[1,1:].sum()=}")
print(f"{total[1]/total[0]=}")

VACC = ('no', 'vacc')


for data_dir in [ '00_lower','01_medium','02_upper']:
    matrices = get_covid_matrices('delta',data_dir,vaccination_statuses=VACC)
    N = matrices['N']
    C =  get_reduced_vaccinated_susceptible_contribution_matrix_covid(R0,'delta',data_dir,vaccination_statuses=VACC)
    yred = get_reduced_eigenvector_covid(R0, 'delta', data_dir,vaccination_statuses=VACC)
    yred = [yred[0], yred[1:].sum()]
    R = C.sum()
    print("======",data_dir,"======")
    print("Reff=",R)
    print(C/R)
    print(C.sum(axis=0))
    print(C.sum(axis=0)/R)
    print("y[uvacc,vacc]=", yred)
    #y = get_reduced_population_eigenvector_covid(R0, 'delta', data_dir)
    y = get_eigenvector_covid(R0,variant='delta',data_dir=data_dir,vaccination_statuses=VACC)

    # ignore young children
    y = y[1:,:]
    # summ vacc status
    newy = np.zeros((3,2))
    newy[:,0] = y[:,0]
    newy[:,1] = y[:,1:].sum(axis=1)

    print(f"{newy/newy.sum(axis=1)[:,None]=}")
    print(f"{newy[1:,:].sum(axis=0)/newy[1:,:].sum()=}")

    print()
    print(f"{newy/newy.sum()=}")
    print(f"{datay=}")
    print()


    Chom =  get_homogeneous_contribution_matrix_covid(R0,'delta',data_dir,vaccination_statuses=VACC)
    Rhom = Chom.sum()
    print(" hom")
    print(Chom/Rhom)
    print(Chom.sum(axis=0))
    print(Chom.sum(axis=0)/Rhom)

    pres = get_circular_vector_and_matrix_presentation(yred, C)
    ax = pres.plot(figsize=(6,6))
    ax.axis('off')
    pres.add_arrows_to_plot()
    pres.add_text_to_plot()

    caps = CircleCapSegmentPresentation(C)
    caps.compute()
    ax = caps.plot()
    ax.get_figure().savefig(f'caps_{data_dir}.pdf')

    C_full = get_contribution_matrix_covid(R0,'delta',data_dir,vaccination_statuses=VACC)

    print(to_ltx(convert_4d_matrix_to_2d_block(C_full/C_full.sum()),
                 fmt='{:4.2f}',
                 separate_columns = range(2,26,2),
                 separate_rows = range(2,26,2),
                 )
          )

#data_dir = '03_realistic'
#
#C = get_reduced_vaccinated_susceptible_contribution_matrix_covid(R0,'delta',data_dir)
#y = get_reduced_eigenvector_covid(R0, 'delta', data_dir)
#y = [y[0], y[1:].sum()]
#R = C.sum()
#print("======",data_dir,"======")
#print(C/R)
#print(C.sum(axis=0))
#print(C.sum(axis=0)/R)
#print(f"{y=}")
#
#

pl.show()

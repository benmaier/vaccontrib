from vaccontrib.covid import (
        get_reduced_vaccinated_susceptible_contribution_matrix_covid,
        get_homogeneous_contribution_matrix_covid,
        get_reduced_eigenvector_covid,
        get_reduced_population_eigenvector_covid,
        get_covid_matrices,
        get_eigenvector_covid,
    )

from vaccontrib.illustration import (
        CircleCapSegmentPresentation,
        get_circular_vector_and_matrix_presentation,
    )

import numpy as np
import matplotlib.pyplot as pl

R0 = [1,1,1,1,1,]
R0 = [1.5,3,3,3,3,]

total = np.array([
            [ 0, 23_119, 189_698, 37_735 ],
            [ 0,  1_102,  66_396, 22_973 ],
        ])

print(total[1,:].sum()/total[0,:].sum())


#for data_dir in ['00_lower', '01_middle', '02_upper','03_realistic']:
for data_dir in ['00_lower', '03_realistic']:
    matrices = get_covid_matrices('delta',data_dir)
    N = matrices['N']
    C =  get_reduced_vaccinated_susceptible_contribution_matrix_covid(R0,'delta',data_dir)
    yred = get_reduced_eigenvector_covid(R0, 'delta', data_dir)
    yred = [yred[0], yred[1:].sum()]
    R = C.sum()
    print("======",data_dir,"======")
    print("Reff=",R)
    print(C/R)
    print(C.sum(axis=0))
    print(C.sum(axis=0)/R)
    print("y[uvacc,vacc]=", yred)
    #y = get_reduced_population_eigenvector_covid(R0, 'delta', data_dir)
    y = get_eigenvector_covid(R0,variant='delta',data_dir=data_dir)

    # ignore young children
    y = y[1:,:]
    # summ vacc status
    newy = np.zeros((3,2))
    newy[:,0] = y[:,0]
    newy[:,1] = y[:,1:].sum(axis=1)

    print(f"{newy/newy.sum(axis=1)[:,None]=}")


    Chom =  get_homogeneous_contribution_matrix_covid(R0,'delta',data_dir)
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

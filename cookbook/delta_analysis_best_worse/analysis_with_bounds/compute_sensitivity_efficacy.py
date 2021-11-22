import numpy as np
import matplotlib.pyplot as pl


from vaccontrib.covid import (
        get_covid_matrices
    )

from vaccontrib.main import (
        get_reduced_vaccinated_susceptible_contribution_matrix,
        get_reduced_vaccinated_susceptible_eigenvector,
        get_eigenvector,
        get_next_generation_matrix_from_matrices,
    )


from tqdm import tqdm

import matplotlib.ticker as mtick

import bfmplot as bp


matrices = get_covid_matrices('delta','01_upper',('no','vacc'))

reduction = np.linspace(1,0,41)

s0 = np.array(matrices['s'])
r0 = np.array(matrices['r'])

n = len(reduction)
Cs = np.zeros((n,2,2))
ys = np.zeros((n,4,2))
breakthroughs = np.zeros((n,3))

#_v = np.array([0.72,0.72,0.72,0.72])
#_v = np.array([0.92,0.92,0.72,0.72])
_v = np.array([1.,1.,1.,1])


for ired, red in enumerate(reduction):
    s = s0.copy()
    r = r0.copy()
    s[:,1] = 1 - (1-s0[:,0] ) * (1-(1-red)*_v)
    r[:,1] = (1-red)*r0[:,0]
    matrices['s'] = s
    matrices['r'] = r

    K = get_next_generation_matrix_from_matrices(1,**matrices)
    C = get_reduced_vaccinated_susceptible_contribution_matrix(K)
    C /= C.sum()
    y = get_eigenvector(K)
    brk = y[1:,:]
    brk = brk[:,1]/brk.sum(axis=1)

    ys[ired,:,:] = y
    Cs[ired,:,:] = C
    breakthroughs[ired,:] = brk

age_groups = ('children','adolescents','adults','elderly')

datavacc = np.array([
            [ 23_119, 159_647, 37_735 ],
            [ 1_102,  66_396, 22_973 ],
        ],dtype=float)

datay = datavacc.T
datay[:,0] -= datay[:,1]
datay = datay / datay.sum()

incidence_distribution_ages = [0.14710685, 0.09388042, 0.60214432, 0.15686842]

x = 1-reduction

fig, ax = pl.subplots(2,2,figsize=(8,8))
ax = ax.flatten()

for u in range(2):
    for v in range(2):
        ax[0].plot(x,Cs[:,u,v],label='uv'[u]+' $\\rightarrow$ '+'uv'[v])

ax[0].legend()
ax[0].set_ylabel('fraction of new infections caused by ...')

data_breakthrough_fractions = datay[:,1]/datay.sum(axis=1)

for i in range(3):
    p, = ax[1].plot(x, breakthroughs[:,i],label=age_groups[i+1])
    ax[1].plot([0.72],[data_breakthrough_fractions[i]],'s',color=p.get_color())

print(datay)
ax[1].set_ylabel('fraction of breakthrough infections in age group')
ax[1].legend()

for i in range(4):
    p, = ax[3].plot(x, ys[:,i].sum(axis=1),label=age_groups[i])
    ax[3].plot([0.72],[incidence_distribution_ages[i]],'s',color=bp.colors[i])

ax[3].legend()

ax[2].plot(x,Cs.sum(axis=1)[:,0],label='unvaccinated')
ax[2].plot(x,Cs.sum(axis=1)[:,1],label='vaccinated')
ax[2].set_ylabel('fraction of new infections caused by ...')
ax[2].legend()
#ax[3].plot(reduction,Cs[:,0,0]+Cs[:,0,1]+Cs[:,1,0])

ax[3].set_ylabel('fraction of new infections within age group')

for a in ax:
    a.set_xlim(0,1)
    a.set_ylim(0,1)
    a.set_xlabel('vaccine efficacy against infection')
    a.yaxis.set_major_formatter(mtick.PercentFormatter(1))
    a.xaxis.set_major_formatter(mtick.PercentFormatter(1))



fig.tight_layout()

fig.savefig('efficacy_curves.pdf')


pl.show()

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
        get_contribution_matrix,
    )


from tqdm import tqdm

import matplotlib.ticker as mtick

import bfmplot as bp

colors = [
            ['#E75740', '#58BDB2'],
            ['#F2957D', '#268D7C'],
        ]

uv_colors = [ colors[0][0], colors[1][1] ]


reduction = np.linspace(1,0,41)
n = len(reduction)

matrices = get_covid_matrices('delta','01_upper',('no','vacc'))
s0 = np.array(matrices['s'])
r0 = np.array(matrices['r'])
b0 = np.array(matrices['b'])

Cs = np.zeros((2,n,2,2))
for imode, reduce_susc_only in enumerate([True,False]):

    _v = np.array([1.,1.,1.,1])

    for ired, red in enumerate(reduction):
        s = s0.copy()
        r = r0.copy()
        b = b0.copy()
        s[:,1] = 1 - (1-s0[:,0] ) * (1-(1-red)*_v)
        if reduce_susc_only:
            r = r0
            b = b0
        else:
            r[:,1] = (1-red)*r0[:,1]
            b[:,1] = (1-red)*b0[:,1] + red * (b0[:,0])
        matrices['s'] = s
        matrices['r'] = r
        matrices['b'] = b

        K = get_next_generation_matrix_from_matrices(1,**matrices)
        C = get_reduced_vaccinated_susceptible_contribution_matrix(K)
        C /= C.sum()
        Cs[imode,ired,:,:] = C

fig, ax = pl.subplots(1,1,figsize=(5,3.5))

x = 1 - reduction

linestyles = ['-','--']

labels = ['const. breakthrough\ntransmissibility reduction',
          'decreasing breakthrough\ntransmissibility reduction',
         ]

ax.plot(x,0.5*np.ones_like(x),c='#aaaaaa',ls='-')
ax.plot([0.22,0.22],[0,.5],c='#aaaaaa',ls='-')
ax.plot([0.41,0.41],[0,.5],c='#aaaaaa',ls='-')

for imode in range(2):
    unvacc = Cs[imode,:,:,:].sum(axis=1)[:,0]
    vacc = Cs[imode,:,:,:].sum(axis=1)[:,1]
    ax.plot(x,unvacc,color=uv_colors[0],label=labels[imode],ls=linestyles[imode])
    ax.plot(x,vacc,color=uv_colors[1],ls=linestyles[imode])
    ax.set_ylabel('fraction of new infections caused by ...')
    ax.legend()
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1))

ax.set_yticks([0,.25,.5,.75,1])
ax.set_xlim(0,1)
ax.set_ylim(0,1)

ax.text(0.85,0.65,'unvaccinated',ha='right',va='top',color=uv_colors[0])
ax.text(0.8,0.1,'vaccinated',ha='right',va='bottom',color=uv_colors[1])

fig.tight_layout()

ax.set_xlabel('age-independent vaccine efficacy s')

bp.strip_axis(ax)


fig.tight_layout()

fig.savefig('efficacy_scan.pdf')


pl.show()

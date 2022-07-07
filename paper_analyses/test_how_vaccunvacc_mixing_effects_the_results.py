import numpy as np
import matplotlib.pyplot as pl


from vaccontrib.covid import (
        get_covid_matrices,
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

colors = [
            ['#E75740', '#58BDB2'],
            ['#F2957D', '#268D7C'],
        ]

uv_colors = [ colors[0][0], colors[1][1] ]


matrices = get_covid_matrices('delta','01_medium,('no','vacc'))

ms = np.linspace(1,0,41)
n = len(ms)

Cs = np.zeros((n, 2,2))
Cunnormeds = np.zeros((n, 2,2))

K = get_next_generation_matrix_from_matrices(1,**matrices)
C = get_reduced_vaccinated_susceptible_contribution_matrix(K)
R0 = 1.2
Rgauge = R0 / C.sum()

for im, m in enumerate(ms):
    mixing = np.array([
                       [ 1.0, m],
                       [ m, 1.0],
                     ])

    K = get_next_generation_matrix_from_matrices(Rgauge,**matrices)
    C = get_reduced_vaccinated_susceptible_contribution_matrix(K)

    K = K * mixing[None,None,:,:]

    C = get_reduced_vaccinated_susceptible_contribution_matrix(K)
    Cunnormeds[im,:,:] = C
    C /= C.sum()
    y = get_eigenvector(K)
    brk = y[1:,:]
    brk = brk[:,1]/brk.sum(axis=1)

    Cs[im,:,:] = C
    print(C)
    print(y)

fig, axs = pl.subplots(1,2,figsize=(8,3))
axs[0].plot(ms, Cs.sum(axis=1)[:,0],label='unvaccinated',c=uv_colors[0])
axs[0].plot(ms, Cs.sum(axis=1)[:,1],label='vaccinated',c=uv_colors[1])
axs[0].legend()
axs[1].plot(ms, Cunnormeds.sum(axis=1)[:,0],label='unvaccinated',c=uv_colors[0])
axs[1].plot(ms, Cunnormeds.sum(axis=1)[:,1],label='vaccinated',c=uv_colors[1])
axs[1].plot(ms, Cunnormeds.sum(axis=1).sum(axis=1),label='total',c='k')
axs[1].legend()


for ax in axs:
    ax.set_xlabel('mixing unvaccinated-vaccinated $m$')
    #ax.set_xlim([0,1])
    #ax.set_ylim([0,1])
    bp.strip_axis(ax)

axs[0].yaxis.set_major_formatter(mtick.PercentFormatter(1))
axs[0].xaxis.set_major_formatter(mtick.PercentFormatter(1))
axs[1].xaxis.set_major_formatter(mtick.PercentFormatter(1))

axs[0].set_ylabel('relative contribution by')
axs[1].set_ylabel('absolute contribution to R by')
#axs[1].set_ylim([0,.5])

fig.tight_layout()
fig.savefig('mixing_analysis.png',dpi=300)

pl.show()




import numpy as np
import matplotlib.pyplot as pl
#pl.rcParams["font.size"] = 12


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

from bfmplot.tools import get_inset_axes

import bfmplot as bp


colors = [
            ['#E75740', '#58BDB2'],
            ['#F2957D', '#268D7C'],
        ]

uv_colors = [ colors[0][0], colors[1][1] ]


matrices = get_covid_matrices('delta','00_lower',('no','vacc'))

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

fig, ax = pl.subplots(1,1,figsize=(3.6,3.3))
ax2 = ax.inset_axes([.7,.3,.3,.2])

axs = [ax2, ax]

x = 1-ms
print(x)

axs[0].plot(x, Cs.sum(axis=1)[:,0],label='unvaccinated',c=uv_colors[0],ls='-.')
axs[0].plot(x, Cs.sum(axis=1)[:,1],label='vaccinated',c=uv_colors[1],ls='--')
axs[1].plot(x, Cunnormeds.sum(axis=1)[:,0],label='unvaccinated',c=uv_colors[0],ls='-.')
axs[1].plot(x, Cunnormeds.sum(axis=1)[:,1],label='vaccinated',c=uv_colors[1],ls='--')
axs[1].plot(x, Cunnormeds.sum(axis=1).sum(axis=1),label='total',c='k',lw=2)
leg = axs[1].legend()
bp.align_legend_right(leg)

axs[1].set_xlim([0,1])
axs[0].set_ylim([0,1])
axs[0].set_xlim([0,1])
axs[1].set_ylim([0,1.2])
axs[1].set_xlabel('mixing unvaccinated-vaccinated m')
axs[0].set_xlabel('m')
axs[1].set_ylabel('absolute contribution to R by')

axs[0].yaxis.set_major_formatter(mtick.PercentFormatter(1))
axs[0].xaxis.set_major_formatter(mtick.PercentFormatter(1))
axs[1].xaxis.set_major_formatter(mtick.PercentFormatter(1))

fig.tight_layout()
axs[0].set_ylabel('relative')
for ax in axs:
    bp.strip_axis(ax)
    ax.set_xticklabels(ax.get_xticklabels()[::-1])

axs[0].set_xticks([0.,1])
axs[0].set_xticklabels(['100%','0%'])
#fig.tight_layout()
#axs[1].set_ylim([0,.5])

fig.savefig('mixing_analysis.png',dpi=150)
fig.savefig('mixing_analysis.pdf',dpi=150)

pl.show()




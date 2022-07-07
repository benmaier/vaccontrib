import numpy as np
import matplotlib.pyplot as pl

from scipy.optimize import root


from compute import make_analysis

from tqdm import tqdm

import matplotlib.ticker as mtick


_1_minus_Ru = np.linspace(0,1,21)

Rv_crit = np.zeros((len(_1_minus_Ru),3))


pl.figure()

def get_root(Ru,dir='01_medium'):
    fun = lambda x: make_analysis([dir],Ru,x,verbose=False)[0]-1
    sol = root(fun,x0=1)

    return sol.x


fig, ax = pl.subplots(1,1,figsize=(3.6,3.3))

vec = []
for i, _1mRu in tqdm(enumerate(_1_minus_Ru)):
    Rv_crit[i,0] = get_root(1-_1mRu,'00_lower')
    Rv_crit[i,1] = get_root(1-_1mRu,'01_medium')
    Rv_crit[i,2] = get_root(1-_1mRu,'02_upper')

for i in range(3):
    p = np.polyfit(_1_minus_Ru, 1-Rv_crit[:,i], 1)
    print("m =", p[0], "-1/m =", -1/p[0])
    ys = np.polyval(p,[0,1])
    vec.append(
        np.array([
            1,
            np.diff(ys),
       ])
    )

for i in range(3):
    vec[i] = np.array([[0,-1],[1,0]]).dot(vec[i])
    vec[i] /= vec[i].sum()

    print(['00_lower','01_medium','02_upper'][i], vec[i])
    print(['00_lower','01_medium','02_upper'][i], vec[i][0]/vec[i][1])

ax.plot(_1_minus_Ru,1-Rv_crit[:,0],color='k',ls=':',label='low')
ax.fill_between(_1_minus_Ru,np.zeros_like(Rv_crit[:,0]), 1-Rv_crit[:,0],color='k',alpha=0.1)
ax.plot(_1_minus_Ru,1-Rv_crit[:,1],color='k',ls='--',label='medium')
ax.fill_between(_1_minus_Ru,np.zeros_like(Rv_crit[:,1]), 1-Rv_crit[:,1],color='k',alpha=0.1)
ax.plot(_1_minus_Ru,1-Rv_crit[:,2],color='k',ls='-',label='high')
ax.fill_between(_1_minus_Ru,np.zeros_like(Rv_crit[:,2]), 1-Rv_crit[:,2],color='k',alpha=0.1)
ax.set_xlim([0,.5])
ax.set_ylim([0,.5])

ax.text(0.03,0.03,'exponential\ngrowth',transform=ax.transAxes,va='bottom',ha='left')
ax.text(0.7,0.6,'epidemic\ncontrol',transform=ax.transAxes,va='bottom',ha='center')

ax.set_xlabel('NPI transmissibility reduction\nunvaccinated')
ax.set_ylabel('NPI transmissibility reduction\nvaccinated')
ax.legend(loc='lower right')

ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
ax.xaxis.set_major_formatter(mtick.PercentFormatter(1))

fig.tight_layout()

fig.savefig('critical_reduction.pdf')
fig.savefig('critical_reduction.png',dpi=300)

import csv

with open('critical_reduction.csv','w') as f:
    writer = csv.writer(f)
    header = (
                'unvacc_zeta_relative_reduction_in_transmissibillity',
                'vacc_zeta_low_eff_relative_reduction_in_transmissibillity',
                'vacc_zeta_med_eff_relative_reduction_in_transmissibillity',
                'vacc_zeta_high_eff_relative_reduction_in_transmissibillity',
            )
    writer.writerow(header)
    for row in zip(_1_minus_Ru, 1-Rv_crit[:,0], 1-Rv_crit[:,1], 1-Rv_crit[:,2]):
        writer.writerow(["{:6.4f}".format(val) for val in row])

pl.show()

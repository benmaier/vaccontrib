import numpy as np
import matplotlib.pyplot as pl

from scipy.optimize import root


from compute import make_analysis

from tqdm import tqdm

import matplotlib.ticker as mtick


_1_minus_Ru = np.linspace(0,1,21)

Rv_crit = np.zeros((len(_1_minus_Ru),2))


pl.figure()

def get_root(Ru,dir='00_lower'):
    fun = lambda x: make_analysis([dir],Ru,x,verbose=False)[0]-1
    sol = root(fun,x0=1)

    return sol.x


fig, ax = pl.subplots(1,1,figsize=(3.6,3.3))

vec = []
for i, _1mRu in tqdm(enumerate(_1_minus_Ru)):
    Rv_crit[i,0] = get_root(1-_1mRu,'00_lower')
    Rv_crit[i,1] = get_root(1-_1mRu,'01_upper')

for i in range(2):
    p = np.polyfit(_1_minus_Ru, 1-Rv_crit[:,i], 1)
    print("m =", p[0], "-1/m =", -1/p[0])
    ys = np.polyval(p,[0,1])
    vec.append(
        np.array([
            1,
            np.diff(ys),
       ])
    )

for i in range(2):
    vec[i] = np.array([[0,-1],[1,0]]).dot(vec[i])
    vec[i] /= vec[i].sum()

    print(['00_lower','01_upper'][i], vec[i])
    print(['00_lower','01_upper'][i], vec[i][0]/vec[i][1])

ax.plot(_1_minus_Ru,1-Rv_crit[:,0],color='k',ls='--',label='lower eff.')
ax.fill_between(_1_minus_Ru,np.zeros_like(Rv_crit[:,0]), 1-Rv_crit[:,0],color='k',alpha=0.1)
ax.plot(_1_minus_Ru,1-Rv_crit[:,1],color='k',ls='-',label='upper bnd.')
ax.fill_between(_1_minus_Ru,np.zeros_like(Rv_crit[:,1]), 1-Rv_crit[:,1],color='k',alpha=0.1)
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

pl.show()

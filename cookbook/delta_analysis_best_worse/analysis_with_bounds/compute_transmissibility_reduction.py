import numpy as np
import matplotlib.pyplot as pl


from compute import make_analysis


Rus = np.linspace(1,.5,101)
Rvs = np.linspace(1,.5,3)

R = np.zeros((len(Rvs),len(Rus),2))


pl.figure()
for iv, Rv in enumerate(Rvs):
    for iu, Ru in enumerate(Rus):
        R[iv, iu, :] = make_analysis(['00_lower','01_upper'],Ru,Rv,verbose=False)

    pl.fill_between(Rus, R[iv,:,0], R[iv,:,1])

pl.show()

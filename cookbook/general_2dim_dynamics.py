import numpy as np

def get_contribution_matrix_eigenvalue_and_eigenvector(K):
    rho, y = get_spectral_radius_and_eigenvector(K)
    C = K * y[None,:]
    return C, rho, y

def get_spectral_radius_and_eigenvector(M):

    shape = M.shape
    assert(shape[0]==shape[1])

    rho, y = np.linalg.eig(M)
    rho = np.abs(rho)

    ndx = np.argmax(rho)

    rho = rho[ndx]
    y = np.abs(y[:,ndx])
    y /= sum(y)

    return rho, y

np.set_printoptions(linewidth=100)
_a = np.array

M, V = 2, 2
K0 = np.arange(1,M**2*V**2+1).reshape(M*V, M*V) / (M**2*V**2*1.5)

C0, R0, y0 = get_contribution_matrix_eigenvalue_and_eigenvector(K0)


quarantine_ndx = [1,2]
K1 = K0.copy()
K1[quarantine_ndx,:] = 0
K1[:,quarantine_ndx] = 0
free_ndx = sorted(list(set(list(range(M*V))) - set(list(quarantine_ndx))))

C1, R1, y1 = get_contribution_matrix_eigenvalue_and_eigenvector(K1)

print("K =\n", K0)
print("C =\n", C0)
print("y_full =\n", y0)
print("y_free =\n", y0[free_ndx]/y0[free_ndx].sum())
print("R_eff =", C0.sum(axis=0))
print("R =", R0)
#print("C_normed =\n", C0/C0.sum())
print()
#print("with R =", R1)
print("K =\n", K1)
print("C =\n", C1)
print("y_full =\n", y1)
print("R_eff =", C1.sum(axis=0))
print("R =", R1)
#print("C_normed =\n", C1/C1.sum())
print()


y = np.array([1,1,1,1.])
growth = [1.]
VV_growth = [1.]
ytot = sum(y)
ytots = [ytot]
y_d_eig = []

for K in [K0,K1]:
    _, eig = get_spectral_radius_and_eigenvector(K)

    for g in range(10):
        y_d_eig.append((y/np.linalg.norm(y)).dot(eig/np.linalg.norm(eig)))
        ytot_old = ytot
        yold = y.copy()
        y = K.dot(y)
        ytot = y.sum()
        ytots.append(ytot)
        growth.append(ytot/ytot_old)
        VV_growth.append(y[free_ndx].sum()/yold[free_ndx].sum())

y_d_eig.append((y/np.linalg.norm(y)).dot(eig/np.linalg.norm(eig)))
import matplotlib.pyplot as pl

fig, ax = pl.subplots(2,2,figsize=(8,4))

i = free_ndx[-1]
ax[0,0].plot(growth)
ax[0,0].plot([0,len(growth)-1],[C1.sum()]*2,':')
ax[0,0].plot([0,len(growth)-1],[C0.sum()]*2,':')

#ax[1,0].plot(VV_growth)
#ax[1,0].plot([0,len(VV_growth)-1],[C1[i,i]]*2,':')
#ax[1,0].plot([0,len(VV_growth)-1],[C0[i,i]]*2,':')

ax[0,1].plot(ytots)
ax[0,0].set_yscale('log')
ax[0,1].set_yscale('log')


ax[1,1].plot(y_d_eig)

pl.show()

import numpy as np
from vaccontrib.main import get_homogeneous_contribution_matrix

v = 0.65
s = 0.72
r = 0.4

C = get_homogeneous_contribution_matrix(1, v, s, r)

print(C/C.sum())
print(v*(1-s)*(1-r)/(1-v))


s = 0.6
C = get_homogeneous_contribution_matrix(1, v, s, r)
print(v*(1-s)*(1-r)/(1-v))

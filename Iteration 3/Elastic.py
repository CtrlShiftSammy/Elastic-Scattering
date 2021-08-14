import math
import cmath
import numpy as np

input_file = open("Iteration 3/elastic.in", "r")
en = (float(input_file.readline().strip()) / 27.211)
rmu = float(input_file.readline().strip())
rstart = float(input_file.readline().strip())
rend0 = float(input_file.readline().strip())
spac = float(input_file.readline().strip())
lmin_str, lmax_str, lspc_str = (input_file.readline()).split(",")
lmin = int(lmin_str.strip())
lmax = int(lmax_str.strip())
lspc = int(lspc_str.strip())
qmod_str, zmod_str, xi_str, eta_str, rmod_str = (input_file.readline()).split(",")
qmod = float(qmod_str.strip())
zmod = float(zmod_str.strip())
xi = float(xi_str.strip())
eta = float(eta_str.strip())
rmod = float(rmod_str.strip())
iborn = int(input_file.readline().strip())
input_file.close()
pi = cmath.pi
ci = complex(0, 1)
con = 0.125
veloc = math.sqrt(2 * en / rmu)
rmu2 = 2.0*rmu
p1 = math.sqrt(2 * en * rmu)
zasy = -1 * qmod
etahyp = float(zasy) / p1
lmin1 = lmin + 1
lmx = lmax + 1
nos = 2* int((rend0-rstart)/2.0/spac) +2
space = (rend0-rstart)/ float(nos)
del1 = space*space*rmu2 /3.0
del2 = 2.0*del1
del4 = 4.0*del1
fl = np.empty(shape = (1 + int((lmx - lmin1) / lspc)))
zl1 = np.empty(shape = (1 + int((lmx - lmin1) / lspc)))

for i in range(lmin1, lmx, lspc):
    fl[i] = float(i * (i-1)) / rmu2
    zl1[i] = 10 ** 20


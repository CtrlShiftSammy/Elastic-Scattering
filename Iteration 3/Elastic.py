import math
import cmath
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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
nmod = int(zmod-qmod+1.0)
zasy =-(zmod-float(nmod)+1.0) #+ 1.d-10
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
ul1 = 0.0
vl1 = 0.0
vvl1 = 0.0
yl1 = 0.0
cen = 0.0
ccen = 0.0
den = 0.0
r = 0.0
v = 0.0
rho = rend0 * p1
tol = 1.0e-04
maximum = 50

def pot_garvey(r):
    global nmod, zmod, eta, xi, zasy, rmod
    if r < rmod:
        v = ((nmod - 1) * (1 - (1.0/((eta/xi)*math.exp(xi*r)-(eta/xi)+1.0))) - zmod) / r
    else:
        v = zasy / r
    return v

fl = np.empty(shape = (1 + lmx))
z11 = np.empty(shape = (1 + lmx))

for i in range(lmin1, lmx + 1, lspc):
    fl[i] = float(i * (i-1)) / rmu2
    z11[i] = 10 ** 20


for i in range(lmin1, lmx + 1, lspc):
    r=rstart
    for j in range(1, nos + 1, 2):
       r=r+space
       v = pot_garvey(r)
       v11=(en-v)*del4
       cen=1.0/(r**2)*del4
       r=r+space
       v = pot_garvey(r)
       vv11=(en-v)*del2
       ccen=1.0/(r**2)*del2
       y11=z11[i]
       w11=v11-cen*fl[i]
       den=1.0+con*w11
       u11=w11/den
       den=1.0+y11
       y11=y11/den-u11
       u11=vv11-ccen*fl[i]
       den=1.0+y11
       z11[i]=y11/den-u11

for i in range(lmin1, lmx + 1, lspc):
    u11=vv11-ccen*fl[i]
    z11[i]=(z11[i]+u11/2.0)/space


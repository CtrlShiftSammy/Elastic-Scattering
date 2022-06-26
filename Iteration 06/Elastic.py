import math
import cmath
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy.lib.type_check import imag
input_file = open("Iteration 06/elastic.in", "r")
en = (float(input_file.readline().strip()) / 27.211)
rmu = float(input_file.readline().strip())
rstart = float(input_file.readline().strip())
ren = float(input_file.readline().strip())
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
etahyp2 = etahyp*etahyp
lmin1 = lmin + 1
lmx = lmax + 1
nos = 2* int((ren-rstart)/2.0/spac) +2
space = (ren-rstart)/ float(nos)
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
rho = ren * p1
tol = 1.0e-04
maximum = 50
def pot(r):
    global nmod, zmod, eta, xi, zasy, rmod
    if r < rmod:
        v = ((nmod - 1) * (1 - (1.0/((eta/xi)*math.exp(xi*r)-(eta/xi)+1.0))) - zmod) / r
    else:
        v = zasy / r
    return v
#data_lines = 182
#data_file = open("Iteration 06/density_xe+2.dat", "r")
data_lines = 270001
data_file = open("Iteration 06/density_xe+2_smooth.dat", "r")
data = np.empty(shape = (data_lines, 8))
for i in range(0, data_lines):
    str1, str2 = (data_file.readline()).strip().split("  ")
    data[i,0] = float(str1.strip())
    data[i,1] = float(str2.strip())
data[0, 2] = 4.0 * pi * ((data[1, 0]) ** 3) * data[1, 1] / 3.0
data[data_lines - 1, 3] = 4.0 * pi * ((data[data_lines - 1, 0]) ** 2) * data[data_lines - 1, 1] / 2.0
data[0, 4] = data[0, 2] / data[0, 0]
for i in range(1, data_lines):
    data[i, 2] = data[i-1, 2] + 4 * pi * ((data[i, 0]) ** 3 - (data[i-1, 0]) ** 3) * (data[i, 1] + data[i-1, 1]) / 6.0
    data[data_lines-i-1, 3] = data[data_lines-i, 3] + 4 * pi * ((data[data_lines-i, 0]) ** 2 - (data[data_lines-i-1, 0]) ** 2) * (data[data_lines-i, 1] + data[data_lines-i-1, 1]) / 4.0
for i in range(1, data_lines):
    data[i, 4] = (zmod - data[i, 2]) / data[i, 0] - data[i, 3]
    data[i, 5] = 0.5 * ((en - data[i, 4]) - ((en - data[i, 4]) ** 2 + 4 * pi * data[i, 1]) ** 0.5)
    data[i, 6] = 0 #add pol pot
    data[i, 7] = data[i, 4] + data[i, 5] + data[i, 6]
print(str(data[data_lines - 1, 2]))
#0 is radius
#1 is rho
#2 is enclosed Q
#3 is pot3
#4 is static pot
#5 is exchange pot
#6 is polarization pot
#7 is total pot
data_file.close()
def pot_data(r):
    global data
    flag = False
    for i in range(0, data_lines):
        if r <= data[i, 0] and flag == False:
            v = data[i-1, 4] + (r - data[i-1, 0]) * (data[i, 4] - data[i-1, 4]) / (data[i, 0] - data[i-1, 0])
            flag = True
            return v
def pot_data_lagrangian(r):
    global data
    sum = 0
    for i in range(0, data_lines):
        term = 1
        for j in range(0, data_lines):
            if i != j:
                term = term * (r - data[j, 0]) / (data[i, 0] - data[j, 0])
                #print(i, j, term)
        sum = sum + term * data[i, 4]
    return sum
def pot_data_smooth(r):
    global data
    m1 = 0
    m2 =0
    c = 0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 0] and r <= data[i+1, 0]):
            if i == 0:
                m1 = (data[i+1, 4] - data[i, 4]) / (data[i+1, 0] - data[i, 0])
                m2 = (data[i+2, 4] - data[i, 4]) / (data[i+2, 0] - data[i, 0])
            elif i == data_lines - 2:
                m1 = (data[i+1, 4] - data[i-1, 4]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+1, 4] - data[i, 4]) / (data[i+1, 0] - data[i, 0])
            else:
                m1 = (data[i+1, 4] - data[i-1, 4]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+2, 4] - data[i, 4]) / (data[i+2, 0] - data[i, 0])
            y = (m2 - m1) * r * r / (2 * (data[i+1, 0] - data[i, 0])) + (m1 * data[i+1, 0] - m2 * data[i, 0]) * r / (data[i+1, 0] - data[i, 0]) + data[i, 4] - ((m2 * data[i, 0] * data[i, 0] - m1 * data[i, 0] * data[i, 0]) / 2 + (m1 * data[i, 0] * data[i+1, 0] - m2 * data[i, 0] * data[i, 0])) / (data[i+1, 0] - data[i, 0])
            return y
def pot_data_quadratic(r):
    global data
    m1 = 0
    m2 =0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 0] and r <= data[i+1, 0]):
            if i == 0:
                m1 = (data[i+1, 7] - data[i, 7]) / (data[i+1, 0] - data[i, 0])
                m2 = (data[i+2, 7] - data[i, 7]) / (data[i+2, 0] - data[i, 0])
            elif i == data_lines - 2:
                m1 = (data[i+1, 7] - data[i-1, 7]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+1, 7] - data[i, 7]) / (data[i+1, 0] - data[i, 0])
            else:
                m1 = (data[i+1, 7] - data[i-1, 7]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+2, 7] - data[i, 7]) / (data[i+2, 0] - data[i, 0])
            a = (m2 - m1) / (2 * (data[i+1, 0] - data[i, 0]))
            b = m1 - 2 * a * data[i, 0]
            c = ((data[i, 7] - b * data[i, 0] - a * data[i, 0] * data[i, 0]) * (data[i + 1, 0] - r) + (data[i+1, 7] - b * data[i + 1, 0] - a * data[i + 1, 0] * data[i + 1, 0]) * (r - data[i, 0])) / ((data[i + 1, 0] - r) + (r - data[i, 0]))
            y = a * r * r + b * r + c
            return y
def pot_final(r):
    if r < rmod:
        return pot_data_quadratic(r)
    else:
        return (qmod) / r
def cgamma(z):
    cof = np.array([76.18009173e0, -86.50532033e0, 24.01409822e0, -1.231739516e0, 1.20858003e-3, -5.36382e-6])
    if z.real < 1:
        ctmp = (z + 5.5) ** (z + 0.5) * cmath.exp(-1 * (z + 5.5))
        cser = 1.000000000190015
        for i in range(1, 7):
            cser = cser + cof[i-1] / (i + z)
        G = cser * ctmp * math.sqrt(2 * pi) / z
        return G
    else:
        ctmp = (z + 4.5) ** (z - 0.5) * cmath.exp(-1 * (z + 4.5))
        cser = 1.000000000190015
        for i in range(0, 6):
            cser = cser + cof[i] / (i + z)
        G = cser * ctmp * math.sqrt(2 * pi)
        return G
def cgammaln(z):
    cof = np.array([57.1562356658629235,-59.5979603554754912,14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5])
    cser = 0.999999999999997092
    ctmp = (z+0.5)*cmath.log(z+5.24218750000000000)-(z+5.24218750000000000)
    for i in range(0, 14): cser = cser + cof[i] / (z + 1 + i)
    lnG = ctmp+cmath.log(2.5066282746310005*cser/z)
    return lnG
def cgammalnphase(z):
    cof = np.array([57.1562356658629235,-59.5979603554754912,14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5])
    cser = 0.999999999999997092
    ctmp = (z+0.5)*cmath.log(z+5.24218750000000000)-(z+5.24218750000000000)
    for i in range(0, 14): cser = cser + cof[i] / (z + 1 + i)
    lnG = ctmp+cmath.log(2.5066282746310005*cser/z)
    phase = imag(lnG)
    while phase >= (2 * math.pi):
        phase -= 2 * math.pi
    return phase
def hyper(i, rho, eta, sigmal):
    thetal = rho - eta*math.log(2*rho) - i*pi*0.5 + sigmal
 
    f1k = 1.0
    g1k = 0.0
    f2k = 0.0
    g2k = 1.0 - eta/rho
 
    f1 = f1k
    g1 = g1k
    f2 = f2k
    g2 = g2k
 
    kount = 0
 
    ak = eta * 0.5 / rho
    bk = (i*(i+1) + eta*eta) * 0.5 / rho
 
    kmin = 1
    kmaximum = 10
 
    kount = kount + 1
 
    for k in range(kmin, kmaximum + 1):        
        f1kp1 = ak*f1k - bk*g1k
        g1kp1 = ak*g1k + bk*f1k
        f2kp1 = ak*f2k - bk*g2k - f1kp1/rho
        g2kp1 = ak*g2k + bk*f2k - g1kp1/rho
        f1 = f1 + f1kp1
        g1 = g1 + g1kp1
        f2 = f2 + f2kp1
        g2 = g2 + g2kp1
        ak = (2.0*k+1) * eta / (2.0*k+2.0) / rho
        bk = (i*(i+1) - k*(k+1) + eta*eta) / (2.0*k+2.0) / rho
        f1k = f1kp1
        g1k = g1kp1
        f2k = f2kp1
        g2k = g2kp1
           
    s = math.sin(thetal)
    c = math.cos(thetal)
 
    fhyp = g1 * c + f1 * s
    ghyp = f1 * c - g1 * s
    dfhyp = g2 * c + f2 * s
    dghyp = f2 * c - g2 * s
    
    return fhyp, ghyp, dfhyp, dghyp
def pl(lin, x):
    pl = 1.0 if lin == 0 else x
    plm2 = 1
    plm1 = x
    for i in range(2, lin + 1):
        pl = ( (2*i-1)*x*plm1 - (i-1)*plm2 ) / i
        plm2 = plm1
        plm1 = pl
    return pl
fl  = np.empty(shape = (1 + lmx))
z11  = np.empty(shape = (1 + lmx))
f  = np.empty(shape = (1 + lmx))
g  = np.empty(shape = (1 + lmx))
fp  = np.empty(shape = (1 + lmx))
gp  = np.empty(shape = (1 + lmx))
phase  = np.empty(shape = (1 + lmx)) 
sigma = np.empty(shape = (1 + lmx))
for i in range(lmin1, lmx + 1, lspc):
    fl[i] = float(i * (i-1)) / rmu2
    z11[i] = 10 ** 20

#schrodinger equation integration
potential_data = np.empty(shape = (nos+1))
for i in range(1, nos + 1):
    potential_data[i] = pot_final(i * space)
    #print(i)

for i in range(lmin1, lmx + 1, lspc):
    r=rstart
    k = 0
    for j in range(1, nos + 1, 2):
       r=r+space
       #v = pot_data_lagrangian(r)
       k+= 1       
       v = potential_data[k]
       v11=(en-v)*del4
       cen=1.0/(r**2)*del4
       r=r+space
       #v = pot_data_lagrangian(r)
       k+= 1
       v = potential_data[k]    
       #print(i, j, v)   
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

for i in range(0, 2):
    cz = i + 1.0 + ci*etahyp
    sigmal = cgammalnphase(cz) 
    f[i], g[i], fp[i], gp[i] = hyper(i, rho, etahyp, sigmal)
for i in range(1, lmax):
    lm1 = i - 1
    lp1 = i + 1
    term1 = ( i + lp1 ) * ( etahyp + ( i * lp1 / rho ) )
    term2 = lp1 * math.sqrt( i * i + etahyp2 )
    term3 = 1.0 / ( i * math.sqrt( lp1 * lp1 + etahyp2 ) )
    f[lp1] = ( term1 * f[i] - term2 * f[lm1] ) * term3
    g[lp1] = ( term1 * g[i] - term2 * g[lm1] ) * term3
    term1 = math.sqrt( lp1*lp1 + etahyp2 ) / lp1
    term2 = ( lp1*lp1 / rho + etahyp) / lp1
    fp[lp1] = term1 * f[i] - term2 * f[lp1]
    gp[lp1] = term1 * g[i] - term2 * g[lp1]
hcon = 4.0*pi/p1/p1
tcs=0.0
phase_file = open("Iteration 06/Phase shifts.txt", "w")
for i in range(lmin1, lmx + 1, lspc):
    l = i - 1
    logder = z11[i] / p1
    fhyp = f[l]
    ghyp = g[l]
    dfhyp = fp[l]
    dghyp = gp[l]
    phase[i] = math.atan( (logder*fhyp-dfhyp)/(dghyp-ghyp*logder) )
    cg = complex(1.0,0.0) 		
    cz = l + 1.0 + ci*etahyp
    sigma[i] = cgammalnphase(cz)
    phase_file.write(str(l) + "   " + str(phase[i]) + "   " + str(sigma[i]) + "\n")
    tcs = tcs + hcon * (2.0*l+1) * math.sin(phase[i])**2    
phase_file.close()
crosssection_file = open("Iteration 06/Cross sections smooth.txt", "w")
for theta in range(0, 180):
    ang = (theta + 0.001) * pi / 180.0
    sint22 = (math.sin(ang / 2)) ** 2
    ruther = -zasy / ( 2.0 * p1*p1 * sint22 )
    ruthzm =  zmod / ( 2.0 * p1*p1 * sint22 )
    cratio = cmath.exp(cgammaln( 1.0+ci*etahyp ) - cgammaln( 1.0-ci*etahyp ))
    cdexpon = cmath.exp( -1*ci*etahyp*math.log(sint22) )
    cfc = ruther * cratio * cdexpon
    cost = math.cos(ang)
    cfnc = complex(0.0,0.0)
    for i in range(lmin1, lmax + 1, lspc):
        l = i - 1
        cfnc = cfnc + (2*l+1.0) * cmath.exp(2.0*ci*sigma[i]) * (cmath.exp(2.0*ci*phase[i])-1.0) * pl(l,cost)
    cfnc = cfnc / (2.0*ci*p1)
    cf = cfc + cfnc
    dcs = abs(cf)**2
    crosssection_file.write(str(theta) + "   " +  str(dcs) + "   " +  str(ruther**2) + "   " + str(ruthzm**2) + "\n")
crosssection_file.close()

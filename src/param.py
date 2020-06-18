#
# model
#
'''
x       : state vector (S,E,I,R,D,N)
dx      : tendency vector (nE,nI,nE2R,nI2R,nD)
p       : parameter vector (sigma,lambd,gamma,delta,mu)
N       : Population
E       : Exposed
I       : Infected
R       : Recovered
D       : Dead
nE      : Newly exposed
nI      : Newly infected
nE2R    : Newly recovered without experiencing any symptoms
nI2R    : Newly recovered after infected
nD      : Newly dead
reproRat: reproduction rate
deathRat: death rate
sigma   : crisis rate (=1/(incubation period); constant)
lambd   : infection rate
gamma   : recovery rate (=1/(infection period); constant)
delta   : mortality rate (=(death rate)/(infection period); constant)
mu      : recovery rate without experiencing any simptoms
'''
reproRat     = [15.,4.5,25.] # min,ave,max
deathRat     = [0.001,0.05,0.20] # min,ave,max
incubaPeriod = [3.,5.,14.] # min,ave,max
#incubaPeriod = [5.,7.,10.] # min,ave,max
infectPeriod = [3.,10.,18.] # min,ave,max
#infectPeriod = [7.,14.,18.] # min,ave,max
sigma  = 1. / incubaPeriod[1]
gamma  = 1. / infectPeriod[1]
delta  = deathRat[1] / infectPeriod[1]
lambd  = reproRat[1] * gamma / 10000000
mu     = 0.01
#mu     = 0
varSigma = 0.05
varGamma = 0.05
varLambd = 100.
varDelta = 0.01
varMu    = 0.
param    = [sigma,lambd,gamma,delta,mu]
varParam = [varSigma,varLambd,varGamma,varDelta,varMu]

#
# observation
#
'''
nobs    : number of observation
yo[0]   : new cases
yo[1]   : new death
obserr  : observation error
'''
nobs = 3
errobs = 10.

#
# data assimilation
#
'''
kmax    : ensemble size
updateParam : whether estimate parameter XXX or not
relaxParam  : whether maintain ensemble spread of parameter XXX or not
prtrbParam  : whether perturb initial value of parameter XXX or not
'''
kmax = 100
undef = -999.
updateParam = [False,True,False,True,False]
relaxParam = [False,True,False,True,False]
prtrbParam = [True,True,True,True,False]

#
# output
#
'''
outNorm : whether output data normalized by population or not
'''
outNorm = True

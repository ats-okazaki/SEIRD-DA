import numpy as np
import param as par
import seird

def LETKF (xin,Y,yo) :
    xa = np.array(xin)
    Xfm, dXf = decompose_x(xin)
    HX, dY = decompose_x(Y)

    #par.nobs = 3 * 10

    # Loop for grid
    for i in range(0,10) :
    #for i in range(0,10) :
        jmax_l = par.nobs
        yo_l  = np.empty(jmax_l,dtype=float)
        dY_l  = np.empty((jmax_l,par.kmax),dtype=float)
        HX_l  = np.empty(jmax_l,dtype=float)
        R_inv = np.zeros((jmax_l,jmax_l),dtype=float)
        Xfm_l = Xfm[i]
        dXf_l = dXf[i,:]
        for j in range(0,par.nobs) :
            yo_l [j]     = yo[j]
            HX_l [j]     = HX[j]
            dY_l [j,:]   = dY[j,:]
            R_inv[j,j] = 1. / (par.errobs * par.errobs)

        # Eigen decomposition
        Pa_tilda = (par.kmax-1.)*np.identity(par.kmax,dtype=float) + np.dot(np.dot(dY_l.T,R_inv),dY_l)
        la, V = np.linalg.eig(Pa_tilda)

        # Mean update
        D = np.zeros((par.kmax,par.kmax),dtype=float)
        for k in range(0,par.kmax) :
            D[k,k] = 1./la[k]   # inv D
        Pa_tilda = np.dot(np.dot(V,D),V.T)
        K = np.dot(np.dot(np.dot(dXf_l,Pa_tilda),dY_l.T),R_inv)
        Xam_l = Xfm_l + np.dot(K,yo_l-HX_l)

        # Ensemble update
        D = np.zeros((par.kmax,par.kmax),dtype=float)
        for k in range(0,par.kmax) :
            D[k,k] = np.sqrt((par.kmax-1.)/la[k])   # D^(-0.5)
        Pa_tilda = np.dot(np.dot(V,D),V.T)
        dXa_l = np.dot(dXf_l,Pa_tilda)

        #
        xa[i,:] = Xam_l + dXa_l[:]

    # negative tendency is not allowed
    for i in range(5) :
        xa[i,:] = np.where( xa[i,:] < 0., 0., xa[i,:])
    # not update params
    for i in range(5) :
        if not par.updateParam[i] :
            xa[5+i,:] = xin[5+i,:]

    return (xa)

#=======================================================
def relax(sprd_f,pa) :
    pout = np.array(pa)
    sprd_a = np.std(pa,ddof=1,axis=1)
    mean_a = np.mean(pa,axis=1)
    for i in range(5) :
        if par.relaxParam[i] :
            if sprd_a[i] != 0. :
                pout[i,:] = mean_a[i] + (pa[i,:] - mean_a[i]) * sprd_f[i] / sprd_a[i]
            else :
                rn = np.random.randn(par.kmax)
                rn = (rn - np.mean(rn)) / np.std(rn,ddof=1)
                #print(rn)
                #print(np.mean(rn))
                #print(np.std(rn,ddof=1))
                pout[i,:] = mean_a[i] + sprd_f[i] * rn[:]
        else :
            pout[i,:] = pa[i,:]
    return (pout)
#=======================================================
def obsope (x1,x0) :
    hx = np.empty((par.nobs,par.kmax),dtype=float)
    # pickup I
    hx[0,:] = x1[2,:] - x0[2,:]
    # pickup D
    hx[1,:] = x1[4,:] - x0[4,:]
    return (hx)
#=======================================================
def decompose_x(xin) :
    xanom = np.array (xin)
    xmean = np.mean (xin, axis=1)
    for k in range(0,par.kmax) :
        xanom[:,k] = xin[:,k] - xmean[:] 
    return (xmean,xanom)
#-------------------------------------------------------
def combine_x(xmean,xanom) :
    xout = np.array(xanom)
    for k in range(0,par.kmax) :
        xout[:,k] = xmean[:] + xanom[:,k]
    return (xout)

#=======================================================
def init (N) :
    xf = np.zeros(6,dtype=float)
    pa = np.zeros(4,dtype=float)
    xf[0] = N
    xf[5] = N
    pa[0] = par.sigma
    pa[2] = par.gamma
    pa[3] = par.delta
    pa[1] = 2.0 * par.gamma / N
    return (xf,pa)
#-------------------------------------------------------
def ensinit (N) :
    xf = np.zeros((6,par.kmax),dtype=float)
    xa = np.zeros((6,par.kmax),dtype=float)
    pa = np.zeros((5,par.kmax),dtype=float)
    # state
    xf[0,:] = N
    xf[5,:] = N
    # sigma
    if par.prtrbParam[0] :
        for k in range(par.kmax) :
            #pa[0,k] = 1. / ( np.random.rand() * (par.incubaPeriod[2]-par.incubaPeriod[0]) + par.incubaPeriod[0] )
            pa[0,k] = 1. / np.random.triangular(par.incubaPeriod[0],par.incubaPeriod[1],par.incubaPeriod[2])
    else :
        pa[0,:] = par.param[0]
    # gamma
    if par.prtrbParam[2] :
        for k in range(par.kmax) :
            #pa[2,k] = 1. / ( np.random.rand() * (par.infectPeriod[2]-par.infectPeriod[0]) + par.infectPeriod[0] )
            pa[2,k] = 1. / np.random.triangular(par.infectPeriod[0],par.infectPeriod[1],par.infectPeriod[2])
    else :
        pa[2,:] = par.param[2]
    # delta
    if par.prtrbParam[3] :
        for k in range(par.kmax) :
            pa[3,k] = np.random.triangular(par.deathRat[0],par.deathRat[1],par.deathRat[2]) / np.random.triangular(par.infectPeriod[0],par.infectPeriod[1],par.infectPeriod[2])
            #while True :
            #    rn = np.random.randn() * par.varParam[3]
            #    if par.param[3] + rn >= 0. :
            #        break
            #pa[3,k] = par.param[3] + rn
    else :
        pa[3,:] = par.param[3]
    # lambda
    if par.prtrbParam[1] :
        for k in range(par.kmax) :
            pa[1,k] = ( np.random.rand() * (par.reproRat[2]-par.reproRat[0]) + par.reproRat[0] ) * par.gamma / N
    else :
        pa[1,:] = par.reproRat[1] * par.gamma / N
    # mu
    #if par.prtrbParam[4] :
    #else :
    pa[4,:] = par.param[4]
    # transform
    pat = transParam( pa )
    return (xf,xa,pa,pat)
#-------------------------------------------------------
def arctransParam(pat) :
    pa = np.array(pat)
    pa[1,:] = np.exp( pat[1,:] )
    pa[3,:] = np.tanh( pat[3,:] * np.pi - 1. ) * 0.5 + 0.5
    return ( pa )
#-------------------------------------------------------
def transParam( pa ) :
    pat = np.array(pa)
    pat[1,:] = np.log(pa[1,:])
    pat[3,:] =  ( np.arctanh((pa[3,:] - 0.5 ) * 2) + 1. ) / np.pi
    return ( pat )
#-------------------------------------------------------
def ensfcst (xin,pa) :
    xout = np.array(xin)
    dx = np.empty((5,par.kmax),dtype=float)
    for k in range(par.kmax) :
        xout[:,k],dx[:,k] = seird.fcst(xin[:,k],pa[:,k])
    return (xout,dx)
#-------------------------------------------------------
def ensupdate (xin,dx) :
    xout = np.array(xin)
    for k in range(par.kmax) :
        xout[:,k] = seird.update(xin[:,k],dx[:,k])
    return (xout)

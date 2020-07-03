import numpy as np

def fcst (xin,pin) :
    #
    # assumption
    #
    '''
    - never infected twice
    - no pre-symptomatic transmission
    '''
    #
    # parameters
    # - p = [sigma,lamda,gamma,delta,mu].T
    #
    '''
    sigma   : crisis rate (1/sigma = incubation periiod; constant)
    lambd   : infection rate
    gamma   : recovery rate (1/gamma = infection period; constant)
    delta   : fatality rate (constant)
    mu      : recovery rate without experiencing symptom
    '''
    #lambd  = 2.0 * gamma / N
    #
    # variables 
    # - x = [S,E,I,R,D,N].T
    # - dx = [nE,nI,nE2R,nI2R,nD].T
    #
    '''
    S : Susceptible
    E : Exposed
    I : Infected
    R : Recoverd
    D : Dead
    N : Population
    '''
    Sin = xin[0]
    Ein = xin[1]
    Iin = xin[2]
    Rin = xin[3]
    Din = xin[4]
    Nin = xin[5]
    sigma = pin[0]
    lambd = pin[1]
    gamma = pin[2]
    delta = pin[3]
    mu    = pin[4]
    #
    # time integration
    #
    nE = min(lambd * Iin * Sin, Sin)
    nI = sigma * Ein
    nE2R = mu * Ein
    nI2R = gamma * Iin
    nD = delta * Iin
    Sout = Sin - nE
    Eout = Ein + nE - nI - nE2R
    Iout = Iin + nI - nI2R -nD
    Rout = Rin + nE2R + nI2R
    Dout = Din + nD
    Nout = Sout + Eout + Iout + Rout + Dout
    #
    # update state but params
    #
    xout = np.array(xin)
    dx = np.empty(5,dtype=float)
    xout[0] = Sout
    xout[1] = Eout
    xout[2] = Iout
    xout[3] = Rout
    xout[4] = Dout
    xout[5] = Nout
    dx  [0] = nE
    dx  [1] = nI
    dx  [2] = nE2R
    dx  [3] = nI2R
    dx  [4] = nD

    return (xout,dx)
#--------------------------------------
def update(xin,dx) :
    Sin = xin[0]
    Ein = xin[1]
    Iin = xin[2]
    Rin = xin[3]
    Din = xin[4]
    Nin = xin[5]
    nE   = dx[0]
    nI   = dx[1]
    nE2R = dx[2]
    nI2R = dx[3]
    nD   = dx[4]
    nE   = max(min(nE, Sin), 0)
    nI   = max(min(nI, Ein), 0)
    nE2R = max(min(nE2R, Ein), 0)
    nI2R = max(min(nI2R, Iin), 0)
    nD   = max(min(nD, Iin), 0)
    Sout = Sin - nE
    Eout = Ein + nE - nI - nE2R
    Iout = Iin + nI - nI2R - nD
    Rout = Rin + nE2R + nI2R
    Dout = Din + nD
    Nout = Sout + Eout + Iout + Rout + Dout
    #
    # update 
    #
    xout = np.empty(np.shape(xin))
    xout[0] = max(Sout,0.)
    xout[1] = max(Eout,0.)
    xout[2] = max(Iout,0.)
    xout[3] = max(Rout,0.)
    xout[4] = max(Dout,0.)
    xout[5] = max(Nout,0.)
    return (xout)

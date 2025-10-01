import numpy as np 
import numba as nb
from numba import njit,prange

def getKinvMat(p,Th,dth,k,k0,k1):
    Axx = np.trapz(p*np.cos(Th)*np.cos(Th),dx=dth,axis=0)
    Ayy = np.trapz(p*np.sin(Th)*np.sin(Th),dx=dth,axis=0)
    Axy = np.trapz(p*np.cos(Th)*np.sin(Th),dx=dth,axis=0)

    kxx = 1/(k+k0*Axx)
    kxy = 1/(k+k0*Axy)
    kyx = 1/(k+k0*Axy)
    kyy = 1/(k+k0*Ayy)

    return (kxx,kxy,kyx,kyy)

@njit
def getNextVelocity_X(nu,nutb,nutc,kb,ubx,ucx,usx,yb,yc,ys,dyb,dyc,dys,kxx,kxy,kyx,kyy,rho,eta,etatc,etatb,phic,phib,uxtop,dt):
    ubxi = 1*ubx
    ucxi = 1*ucx
    usxi = 1*usx

    # fluid domain
    usxi[-1] = uxtop # top boundary condition
    usxi[1:-1] = usxi[1:-1] + dt*nu*(usxi[2:] - 2*usxi[1:-1] + usxi[0:-2] )/(dys*dys)
    usxi[0] = (eta*usx[1] + (etatc/phic)*ucx[-2])/(eta+etatc/phic) # bot boundary condition

    # cartilage domain 
    ucxi[-1] = usxi[0] # top boundary condition
    ucxi[1:-1] = ucxi[1:-1] + dt*nutc*(ucxi[2:] - 2*ucxi[1:-1] + ucxi[0:-2] )/(dyc*dyc) - dt*nu*kxx[1:-1]*ucxi[1:-1]
    ucxi[0] = ((etatc/phic)*ucx[1] + (etatb/phib)*ubx[-2])/(etatc/phic + etatb/phib) # bot boundary condition

    # bone domain 
    ubxi[-1] = ucxi[0] # top boundary condition
    ubxi[1:-1] = ubxi[1:-1] + dt*nutb*(ubxi[2:] - 2*ubxi[1:-1] + ubxi[0:-2] )/(dyb*dyb) - dt*nu*kb*ubxi[1:-1]
    ubxi[0] = 0.0 # bot boundary condition

    return (ubxi, ucxi, usxi)

def interectionFunc(vmag,phi,Th,ths):
    sigdisp = 1/(1+30*vmag)
    Ph = np.tile(phi,(len(ths),1))
    sigdispinv = 1/sigdisp
    return 0.5*np.exp(sigdispinv*np.cos(Th-Ph)) + 0.5*np.exp(sigdispinv*np.cos(Th-Ph+np.pi))

def d2pdth2(p,dth):
    secondDerivative = (np.roll(p,1,axis=0) + np.roll(p,-1,axis=0) - 2*p)/(dth*dth)
    return secondDerivative

def calcMuSig(p,Th,yc):
    cs = p*np.cos(Th)
    sn = p*np.sin(Th)
    
    mu = np.trapz(p*Th,axis=0)/np.pi
    mu2 = np.tile(mu,(len(ths),1))
    sig = np.std(p*(Th-mu2)*(Th-mu2),axis=0)/np.pi

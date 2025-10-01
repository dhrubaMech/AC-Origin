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
def getNextVelocity_X(nu,nutc,nutb,kb,ubx,ucx,usx,yb,yc,ys,dyb,dyc,dys,kxx,kxy,kyx,kyy,rho,eta,etatc,etatb,phic,phib,uxtop,dt):
	ubxi = 1*ubx
	ucxi = 1*ucx
	usxi = 1*usx

	# fluid domain
	usxi[-1] = uxtop # top boundary condition
	usxi[1:-1] = usxi[1:-1] + dt*nu*(usxi[2:] - 2*usxi[1:-1] + usxi[0:-2] )/(dys*dys)
	usxi[0] = (eta*usx[1] + (etatc/phic)*ucx[-2])/(eta+etatc/phic) # bot boundary condition

	# cartilage domain 
	ucxi[-1] = usxi[0] # top boundary condition
	#print(ucxi[1:-1].shape , kxx[1:-1].shape  )
	ucxi[1:-1] = ucxi[1:-1] + dt*nutc*(ucxi[2:] - 2*ucxi[1:-1] + ucxi[0:-2] )/(dyc*dyc) - dt*nu*kxx[1:-1]*ucxi[1:-1]
	ucxi[0] = ((etatc/phic)*ucx[1] + (etatb/phib)*ubx[-2])/(etatc/phic + etatb/phib) # bot boundary condition

	# bone domain 
	ubxi[-1] = ucxi[0] # top boundary condition
	ubxi[1:-1] = ubxi[1:-1] + dt*nutb*(ubxi[2:] - 2*ubxi[1:-1] + ubxi[0:-2] )/(dyb*dyb) - dt*nu*kb*ubxi[1:-1]
	ubxi[0] = 0.0 # bot boundary condition

	return (ubxi, ucxi, usxi)

def interectionFunc(vmag,phi,Th,ths):
	sigdisp = 1/(1+70*vmag)
	Ph = np.tile(phi,(len(ths),1))
	sigdispinv = 1/sigdisp
	return 0.5*np.exp(sigdispinv*np.cos(Th-Ph)) + 0.5*np.exp(sigdispinv*np.cos(Th-Ph+np.pi))

def d2pdth2(p,dth):
	secondDerivative = (np.roll(p,1,axis=0) + np.roll(p,-1,axis=0) - 2*p)/(dth*dth)
	return secondDerivative

def calcMuSig(p,Th,yc,ths):
	cs = p*np.cos(Th)
	sn = p*np.sin(Th)
	
	mu = np.trapz(p*Th,axis=0)/np.pi
	mu2 = np.tile(mu,(len(ths),1))
	sig = np.std(p*(Th-mu2)*(Th-mu2),axis=0)/np.pi


#@njit
#def generateVelocities(ts,Tpswing,Tpnormal,us0,un0,overlap):
#    # https://asmedigitalcollection.asme.org/tribology/article/122/1/332/477319/Contact-Fatigue-Failure-of-Ultra-High-Molecular
#    Tpcycle = Tpnormal+Tpswing
#    tint = ts % Tpcycle
#    if tint <= Tpnormal:
#        # loading
#        ux0 = 0
#        uy0 = un0*np.sin(2*np.pi*tint/(Tpnormal*(1+overlap)))
#    #elif Tpnormal < tint < Tpnormal*(1+overlap):
#    elif tint >= Tpnormal*(1+overlap):
#        # swinging
#        ux0 = us0*np.sin(2*np.pi*tint/(Tpswing))
#        uy0 = 0.0
#    else:
#    	# overlap
#    	ux0 = us0*np.sin(2*np.pi*tint/(Tpswing))
#    	uy0 = un0*np.sin(2*np.pi*tint/(Tpnormal*(1+overlap)))
#    return (ux0,uy0)

@njit
def generateVelocities(ts,Tpswing,Tpnormal,us0,un0):
    # https://asmedigitalcollection.asme.org/tribology/article/122/1/332/477319/Contact-Fatigue-Failure-of-Ultra-High-Molecular
    Tpcycle = Tpnormal+Tpswing
    tint = ts % Tpcycle
    if tint < Tpnormal:
        # loading
        ux0 = 0
        uy0 = un0*np.sin(2*np.pi*tint/(Tpnormal))
    else:
        # swinging
        ux0 = us0*np.sin(2*np.pi*tint/(Tpswing))
        uy0 = 0.0
    return (ux0,uy0)

@njit(parallel=True)
def calcFibDist(NbinH,phfib,phfibdist,bins,binsize,labels):
	for i in prange(NbinH):
		idx2 = np.where(labels==i)[0]
		if idx2.shape[0] > 0:
			hist, bin_edges = np.histogram(phfib[idx2],bins=bins)
			
			#pint = np.trapz(hist,dx=binsize) ; hist = hist/pint		## numerically more accurate/precise
			#hist[np.isnan(hist)] = 0.0
			
			#hist = hist*binsize/sum(hist*binsize)
			
			#bin_widths = bin_edges[1:] - bin_edges[:-1]
			#hist = hist / (np.sum(hist) * bin_widths)
			
			phfibdist[:,i] = hist
	#print(phfibdist) ; exit()
	return phfibdist

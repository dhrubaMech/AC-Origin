import numpy as np 
import matplotlib.pyplot as plt 
import numba as nb
from numba import njit,prange

from param import *

from matplotlib import rc
plt.rcParams["font.family"] = "serif"   
plt.rcParams["mathtext.fontset"] = "stix"
rc('text', usetex=True)
tickFontSize=14
labelFontSize=18
figPanelFontSize=24

from moduleFluidProp import getNextVelocity_X,getKinvMat,interectionFunc,d2pdth2

def simulate():
    from param import ubx, ucx, usx
    
    p = 0*Y + 1.0/np.pi
    (kxx,kxy,kyx,kyy) = getKinvMat(p,Th,dth,k,k0,k1)
    f, ax = plt.subplots(2,2,figsize=(8,8))  
       
    for ti, ts in enumerate(t):
        print(ti)
        (ubx, ucx, usx) = getNextVelocity_X(nu,nutb,nutc,kb,ubx,ucx,usx,yb,yc,ys,dyb,dyc,dys,kxx,kxy,kyx,kyy,rho,eta,etatc,etatb,phic,phib,ux0[ti],dt)
        ucy = -(uy0[ti]/kyy)*(1.1*hc - yc/hc)

        phi = np.arctan2(ucy,ucx)
        phi[phi<-np.pi/2] += np.pi
        phi[phi>np.pi/2] -= np.pi
        
        u = np.sqrt(ucx*ucx+ucy*ucy)

        # fiber evolution
        dpdt = alp*(1+eps*interectionFunc(u,phi,Th,ths)) + D*d2pdth2(p,dth) - bet*p
        dpdt += nois*np.random.normal(loc=1,scale=1,size=p.shape)

        p += dt*dpdt 
        pint = np.trapz(p,dx=dth,axis=0)
        pint = np.tile(pint,(len(ths),1))
        p = p/pint
        # fiber evolution done

        (kxx,kxy,kyx,kyy) = getKinvMat(p,Th,dth,k,k0,k1)
        if ti % 500 == 499:
            phi = np.arctan2(ucy,ucx)
            phi[phi<-np.pi/2] += np.pi
            phi[phi>np.pi/2] -= np.pi

            ax[0,0].plot(ucx,yc,'g')
            ax[0,0].plot(ucy,yc,'g--')

            ax[0,0].plot([0,0],[0,1],'k--',lw=1)
            ax[0,0].contourf(Th,Y,p,100)#np.linspace(0,0.33,100))
            ax[0,0].contourf(Th+np.pi-dth,Y,p,100)#np.linspace(0,0.33,100))
            ax[0,0].set_xlim([-np.pi/2,np.pi/2])

            ax[0,1].plot(kxx,yc,'r')
            ax[0,1].plot(kyy,yc,'g')

            ax[1,0].plot(ths,p[:,0],'r')
            ax[1,0].plot(ths,p[:,9],'g')
            ax[1,0].plot(ths,p[:,14],'b')
            ax[1,0].plot(ths,p[:,19],'k')
            
            ax[1,1].axis("off")
            
            ax[0,0].set_xlabel("$\\theta$",fontsize=labelFontSize)
            ax[0,0].set_ylabel("$y$",fontsize=labelFontSize)
            
            ax[0,1].set_xlabel("$\mathbf{K}_{xx}$ and $\mathbf{K}_{yy}$",fontsize=labelFontSize)
            ax[0,1].set_ylabel("$y$",fontsize=labelFontSize)
            
            ax[1,0].set_xlabel("$\\theta$",fontsize=labelFontSize)
            ax[1,0].set_ylabel("$y$",fontsize=labelFontSize)
            
            for axi in ax.flatten():
            	axi.tick_params(labelsize=tickFontSize)

            plt.pause(0.5)
            ax[0,0].cla()
            ax[0,1].cla()
            ax[1,0].cla()
            ax[1,1].cla()
            
            
    plt.tight_layout()
    plt.show()
    

if __name__ == "__main__":
    simulate()


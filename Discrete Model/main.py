import numpy as np , time, os
import pickle as pkl
from PIL import Image
import pyvista as pv

from paramActivity import *

from modulePlotting import getCMapValues,renderfiber,plotfibers,pltDist
from moduleFluidProp import getKinvMat,generateVelocities,getNextVelocity_X,interectionFunc,calcFibDist
from moduleACgeneration import genCollagenCM

def simulate():
	#from param import ubx, ucx, usx
	
	Tps = [0.05,0.95,0.5,0.05,0.95]		## Extension/Flexion time period for different activities
	Tpn = [0.05,0.05,0.5,0.95,0.95]     ## Normal loading time period for different activities
	
	for stype in range(len(stype)):
		for sm in range(Nsample):
			path = f"Data/Activities/sample{sm}/"
			if not os.path.isdir(path):
				os.makedirs(path)
			
			Tpswing  = Tps[stype]
			Tpnormal = Tpn[stype]
			
			ubx = 0*yb
			ucx = 0*yc
			usx = 0*ys
			
			p = 0*Y + 1.0/np.pi
			(kxx,kxy,kyx,kyy) = getKinvMat(p,Th,dth,k,k0,k1)
			#print(kxx.shape)
			
			xf0 = [] ; yf0 = [] ; zf0 = [] ; ph0fib = [] ; binIDz = [] ; th0fib = []
			
			ncolrs = len(ths)
			phspace = np.linspace(-np.pi/2,np.pi/2,ncolrs+1)
			#phspace = np.linspace(0,np.pi/2,ncolrs+1)
			colname = "hsv" #"hsv"				# or "custom"
			colors = getCMapValues(colname,ncolrs)
			colIdx = [] ; colfibs = [] 
			
			bins,binsize = np.linspace(-np.pi/2,np.pi/2,NbinPhi+1,endpoint=True,retstep=True)
			phfibdist = np.zeros((NbinPhi,NbinH))
			
			" for plottings "
			viewdist = False
			if viewdist:
				f, ax = plt.subplots(2,2,figsize=(6,6))	 
			viewLive = False
			if viewLive:
				fig,ax = plt.subplots(1,4,figsize=(16,4.0))
			viewLast = False
			
			filename = f"Type{stype+1}_AC_alphat{alpha}_taubeta{tau_beta}_us{us0}_un{un0}_Tps{Tpswing}_Tpn{Tpnormal}_Tp{Tp}_dt{dt}"
			print(filename)
			
			begin = time.time()
			kappa_n0 = 0 ; ti = 0 ; ts = dt
			
			pl = pv.Plotter(window_size=[800,800],lighting='none',off_screen=True)
			fibdata = pv.PolyData()
			light = pv.Light()
			light.set_direction_angle(20,-20)
			pl.add_light(light)
			
			for ti, ts in enumerate(t[1:]):
			#while True:
				#print(ti)
				#(ux0,uy0) = generateVelocities(ts,Tpswing,Tpnormal,us0,un0,overlap)
				(ux0,uy0) = generateVelocities(ts,Tpswing,Tpnormal,us0,un0)
				(ubx, ucx, usx) = getNextVelocity_X(nu,nutc,nutb,kb,ubx,ucx,usx,yb,yc,ys,dyb,dyc,dys,kxx,kxy,kyx,kyy,rho,eta,etatc,etatb,phic,phib,ux0,dt)
				ucy = -(uy0/kyy)*(1.1*hc - yc/hc)

				# fiber evolution
				
				phi = np.arctan2(ucy,ucx)
				phi[phi < -np.pi/2] += np.pi
				phi[phi >  np.pi/2] -= np.pi
				u = np.sqrt(ucx*ucx+ucy*ucy)
				thsdist = interectionFunc(u,phi,Th,ths)
				
				pint = np.trapz(thsdist,dx=dth,axis=0)
				pint = np.tile(pint,(len(ths),1))
				thsdist = thsdist/pint

				" fiber generation "
				Ngen = int(np.random.choice(alpdist))
				xf,yf,zf,ph0f,idz = genCollagenCM(Wsup,Lsup,Hsup,Ngen,yc,dth,ths,thsdist)
				thsf = np.random.choice(thetaDist,size=Ngen)
				xf0 = np.concatenate((xf0,xf))
				yf0 = np.concatenate((yf0,yf))
				zf0 = np.concatenate((zf0,zf))
				ph0fib = np.concatenate((ph0fib,ph0f))
				th0fib = np.concatenate((th0fib,thsf))
				binIDz = np.concatenate((binIDz,idz))
				
				" color of the fiber "
				cid = np.zeros(Ngen,dtype=int)
				for i in range(Ngen):
					cid[i] = np.argmin(abs(phspace-ph0f[i]))
				colIdx = np.concatenate((colIdx,cid))
					
				" fiber removal "
				Nfibs_ = len(xf0)
				betai  = np.random.rand(Nfibs_)
				#remove = np.where((1-betai) < tau_beta)[0]
				remove = np.where(betai > tau_beta)[0]
				xf0 = np.delete(xf0,remove)
				yf0 = np.delete(yf0,remove)
				zf0 = np.delete(zf0,remove)
				th0fib = np.delete(th0fib,remove)
				ph0fib = np.delete(ph0fib,remove)
				binIDz = np.delete(binIDz,remove)
				colIdx = np.delete(colIdx,remove)
				Nfibs  = len(xf0)
				
				# fiber evolution done
				
				" calc fiber-orientation distribution "
				phfibdist = calcFibDist(NbinH,ph0fib,phfibdist,bins,binsize,binIDz)
				
				" calc kappa "
				if ti > interval:
					pint = np.trapz(phfibdist,dx=dth,axis=0)
					pint = np.tile(pint,(len(ths),1))
					phfibdist = phfibdist/pint
					(kxx,kxy,kyx,kyy) = getKinvMat(phfibdist,Th,dth,k,k0,k1)
					
					kappa_n1 = (kxx- kyy)/(kxx+kyy)
					deltaKappa = np.max(np.abs(kappa_n1 - kappa_n0))
					kappa_n0 = kappa_n1
					
					if (((ti+1)%interval == 0) or ti==0):
						print(f"Iter {ti+1}/{t.shape[0]} | {np.round((ti/t.shape[0])*100,2)}% | Nfibs {Nfibs} | dKappa {deltaKappa}",remove.shape[0])
						
						renderfiber(pl,fibdata,xf0,yf0,zf0,Nfibs,ph0fib,th0fib,colors,colIdx,meanLc,f"{ti}")
						if viewLive:
							plotfibers(xf0,zf0,ph0fib,Nfibs,colors,colIdx,ti,Wsup,Hsup,phfibdist,kxx,kyy,ax)
						if viewdist:
							pltDist(ucy,ucx,p,kxx,kyy)
						if deltaKappa < 1e-10:
							break
				#ti += 1 ; ts += dt
				exit()
			
			print(f"Total fibers {len(xf0)} | total nodes {len(xf0)*discNp} | time elapsed : {np.round((time.time()-begin)/60,3)} mins")
			
			filename += f"_Nfibs{len(xf0)}_Np{len(xf0)*discNp}"
			pkl.dump((xf0,yf0,zf0,ph0fib,phfibdist,kxx,kxy,kyx,kyy,yc,colors,colIdx),open(f"{path}{filename}.pkl","wb"))
			
			if viewLive:
				plt.savefig(f"{path}{filename}.png",dpi=300)
				plt.close()
				
			if viewLast:
				fig,ax = plt.subplots(1,4,figsize=(16,4.0))
				#plotfibers(xf0,zf0,Nfibs,colors,colIdx,ti,t,Wsup,Hsup,phfib,phfibdist,kxx,kyy,ax)
				plotfibers(xf0,zf0,ph0fib,Nfibs,colors,colIdx,ti,Wsup,Hsup,phfibdist,kxx,kyy,ax)
				plt.savefig(f"{path}{filename}_LastTimeStep.png",dpi=300)
				plt.close()
		
	


if __name__ == "__main__":

	simulate()



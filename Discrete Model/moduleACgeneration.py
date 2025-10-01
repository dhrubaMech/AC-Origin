import sys
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

#np.random.seed(999)

from matplotlib import rc
plt.rcParams["font.family"] = "serif"   
plt.rcParams["mathtext.fontset"] = "stix"
rc('text', usetex=True)
tickFontSize=16
labelFontSize=20
figPanelFontSize=24

from modulePDFs import generatePointsOnSphere,getNormalPDF,getpoints,getFibDistr,getRandomOnCellSurface

@njit
def genCollagenCM(Wsup,Lsup,Hsup,Ngen,zcs,dphs,phs,phipdf):
	"<--- gen seed points --->"
	x0 = Wsup*np.random.rand(Ngen)
	y0 = Lsup*np.random.rand(Ngen)
	z0 = Hsup*np.random.rand(Ngen)

	"<--- th0 of the fiber --->"
	#th0fib = np.random.choice(thetaDist)		## this I can choose later too
	
	"<--- ph0 of the fiber --->"
	ph0fib = np.zeros(Ngen)
	idz = np.zeros(Ngen)
	for i in range(Ngen):
		iz = np.argmin(np.abs(zcs - z0[i])) 
		idz[i] = iz
		pdf = phipdf[:,iz]	
		
		area = np.trapz(pdf,dx=dphs)
		pdf = pdf/area
		#print(np.trapz(pdf,dx=dphs),np.allclose(np.sum(pdf), 1.0))
		#exit()
		
		pdf = pdf/np.sum(pdf)		### forcing the sum to be exactly 1
		#ph0fib[i] = np.random.choice(phs,p=pdf)
		ph0fib[i] = phs[np.searchsorted(np.cumsum(pdf), np.random.random(), side="right")]
	
	"<--- Creating steps/Walks for each angles --->" ## this I can draw during plotting
	#xf = x0 + fibLen*np.cos(ph0fib)
	#yf = y0 + fibLen*np.sin(th0fib)
	#zf = z0 + fibLen*np.sin(ph0fib)    
	
	return x0,y0,z0,ph0fib,idz

@njit
def genCollagenCMOA(Wsup,Lsup,Hsup,Ngen,zcs,dphs,phs,phipdf,zone):
	"<--- gen seed points --->"
	x0 = Wsup*np.random.rand(Ngen)
	y0 = Lsup*np.random.rand(Ngen)
	z0 = Hsup*np.random.rand(Ngen)
	dist = np.sqrt((x0-zone[0])**2 + (z0-zone[1])**2)

	"<--- th0 of the fiber --->"
	#th0fib = np.random.choice(thetaDist)		## this I can choose later too
	
	"<--- ph0 of the fiber --->"
	ph0fib = np.zeros(Ngen)
	idz = np.zeros(Ngen)
	for i in range(Ngen):
		iz = np.argmin(np.abs(zcs - z0[i])) 
		idz[i] = iz						## storing the z bin postion
		
		" noise "
		if dist[i] < zone[2]:
			pdf = np.ones(phipdf[:,iz].shape[0])*0.1
		else:
			pdf = phipdf[:,iz]
#		pdf = phipdf[:,iz]
			
		pdf = pdf*dphs/sum(pdf*dphs)
		#ph0fib[i] = np.random.choice(phs,p=pdf)		### storing the orientation of the fiber
		ph0fib[i] = phs[np.searchsorted(np.cumsum(pdf), np.random.random(), side="right")]
	
	return x0,y0,z0,ph0fib,idz
	









	
	

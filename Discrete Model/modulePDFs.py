import numpy as np, sys
from numba import jit
from numba.typed import List
import matplotlib.pyplot as plt

#@jit(nopython=True)
def generatePointsOnSphere(ph0,th0,N,angWidth):
	us = np.random.rand(N)
	
	vmin = (np.cos(angWidth)+1)/2
	vs = (1 - vmin)*np.random.rand(N) + vmin

	ths = 2*np.pi*us
	phs = np.arccos(2*vs-1)

	xs = np.sin(phs)*np.cos(ths)
	ys = np.sin(phs)*np.sin(ths)
	zs = np.cos(phs)	

	# rotation matrix calculation
	# https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
	a = np.array([0.0,0.0,1.0])
	b = np.array([np.sin(ph0)*np.cos(th0),np.sin(ph0)*np.sin(th0),np.cos(ph0)])
	v = np.cross(a,b)
	vx = np.array([[0, -v[2], v[1]],
					[v[2], 0, -v[0]],
					[-v[1], v[0], 0]])
	R = np.eye(3) + vx + np.dot(vx,vx)/(1+np.dot(a,b))
	points = np.vstack((xs,ys,zs))#.transpose()
	rotatedPoints = np.dot(R,points)
	xs = rotatedPoints[0,:]
	ys = rotatedPoints[1,:]
	zs = rotatedPoints[2,:]
	
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')
	# ax.scatter3D(xs,ys,zs,s=1,c="r")
	# plt.show()
	
	return (xs,ys,zs)


@jit(nopython=True)  
def getNormalPDF(th0,thw):
	# gaussian
	ths = np.linspace(-np.pi,np.pi,500)
	pdf = (1.0/np.sqrt(2*np.pi*thw*thw))*np.exp(-((ths - th0)**2)/(2*thw*thw))
	pdf = pdf/np.trapz(pdf,ths)

	#plt.plot(ths*180/np.pi,pdf,'.-'); plt.ylim([0,1.2*max(pdf)]); plt.show();# exit()   
	dths = ths[1]-ths[0]
	pdf = pdf*dths/sum(pdf*dths)
	
	#plt.plot(ths/np.pi,pdf,'.-'); plt.ylim([0,1.2*max(pdf)]); plt.show();sys.exit()
	return (ths,pdf)

def getpoints(W,L,H,std):
	x = np.random.normal(W/2,scale=std,size=1)
	y = np.random.normal(L/2,scale=std,size=1)
	z = np.random.normal(H/2,scale=std,size=1)
	
	if (0 <= x and x <= W):
		if (0 <= y and y <= L):
			if (0 <= z and z <= H):
				return x,y,z
			else:
				return getpoints(W,L,H,std)
		else:
			return getpoints(W,L,H,std)
	else:
		return getpoints(W,L,H,std)

def getSeedPointsExpo(W,L,H,Nsub=int(1e4),N=int(1e5)):
	x = np.random.rand(N)*W
	y = np.random.rand(N)*L
	z = np.random.rand(N)*H
	
	dist = np.sqrt((x-W/2)**2+(y-L/2)**2+(z-H/2)**2)
	a = 0.1 ; b = 1.0 ; c = 0.5 ;
	pdf = c + np.exp(b - dist/a)
	
	sort = np.argsort(dist)
	plt.plot(dist[sort],pdf[sort])
	plt.show()
	exit()
	
	return 

def getSeedPointsPowerlaw(W,L,H,n,Nsub=int(1e4),N=int(1e5)):
	x = np.random.rand(N)*W
	y = np.random.rand(N)*L
	z = np.random.rand(N)*H
	
	dist	= np.sqrt((x-W/2)**2+(y-L/2)**2+(z-H/2)**2)
	maxdist = np.sqrt((W-W/2)**2+(L-L/2)**2+(H-H/2)**2) 
	pdf = 1 - (dist/maxdist)**-n
	normpdf = pdf/np.sum(pdf)
	
	index = np.random.choice(np.arange(x.shape[0]),size=Nsub,replace=False,p=normpdf)
	x = x[index]
	y = y[index]
	z = z[index]
	
	return (x,y,z) 

def getFibDistr(W,L,H,center,cellrad,fibrad,fibAngleSTD,effLen,distrSlope,Npoints=int(1e7)):
	### effLen 	  = length of decay [length upto which fibers are dense]
	### lowestpdf = pdf saturation [if the value is close to "highpdf" the distribution will be uniform and if the value is very low the distribution will be dense around the cell only]
	### highpdf   = max value of pdf
	effLen = 1*effLen
	
	x = np.random.rand(Npoints)*W ; y = np.random.rand(Npoints)*L ; z = np.random.rand(Npoints)*H
	dist  = np.sqrt((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)
	index = np.where(dist > (cellrad+fibrad))[0]
	x = x[index] ; y = y[index] ; z = z[index] ; dist = dist[index]

	#lowestpdf = 1e-2 ; highpdf = 1.0 ;
	lowestpdf = 0.1 ; highpdf = 1.0 ;  
	fibDensity= sigmoidPDF(dist,effLen+cellrad,distrSlope,1)
	#fibDensity= sigmoidPDF(dist,1.8*cellrad,distrSlope,1)
	fibDensity  = normalizeVals(fibDensity,lowestpdf,highpdf)
	normpdf   = fibDensity/np.sum(fibDensity)
	
	if False:
		fig = plt.figure(figsize=(7,7))
		ax = fig.add_subplot(111, projection='3d')
		seed = np.random.choice(np.arange(x.shape[0]),size=int(1e4),replace=False,p=normpdf)
		ax.scatter3D(x[seed],y[seed],z[seed],s=3,c="k",alpha=0.5)
		ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
		ax.set_box_aspect((1,1,1))
		plt.show()
		exit()
	
	if False:
		sort = np.argsort(dist)
		#plt.plot(dist[sort],normpdf[sort])
		seed = np.random.choice(np.arange(x.shape[0]),size=int(1e4),replace=False,p=normpdf)
		theta = np.linspace(0,2*np.pi)
		plt.plot(center[0]+cellrad*np.cos(theta),center[1]+cellrad*np.sin(theta),lw=1,ls="--",c="r")
		plt.scatter(x[seed],y[seed],s=0.5,c="k",marker="o")
		#plt.scatter(x,y,s=1,c="k")
		plt.axis("equal")
		plt.show()
		exit()
	
	AngleSTD = sigmoidPDF(dist,(effLen+cellrad)*1.2,distrSlope*3,-1)
	#AngleSTD = sigmoidPDF(dist,1.8*cellrad,distrSlope,-1)
	AngleSTD = normalizeVals(AngleSTD,0,np.pi)
	
	if False:
		index = np.argsort(dist) ; 
		fig, ax1 = plt.subplots()
		ax2 = ax1.twinx()
		ax1.scatter(dist[index],fibDensity[index],c="r",s=2)
		#ax2.scatter(dist[index],AngleSTD[index]*180/np.pi,c="k",s=2)
		
		ax1.plot([cellrad+fibrad,cellrad+fibrad],[min(fibDensity),max(fibDensity)],c="r",ls="--")
		ax1.plot([cellrad,cellrad],[min(fibDensity),max(fibDensity)],c="g",ls="--")
		ax1.plot([cellrad+effLen,cellrad+effLen],[min(fibDensity),max(fibDensity)],c="b",ls="--")
		
		ax1.set_xlim([0,3*cellrad])
		ax1.set_ylabel("p",c="r")
		ax2.set_ylabel("$\\theta$",c="k")
		plt.tight_layout()
		plt.show()
		exit()
	
	index = np.where(dist < effLen+cellrad)[0]
	x = x[index] ; y = y[index] ; z = z[index] ; 
	normpdf = normpdf[index]/np.sum(normpdf[index]) ; AngleSTD = AngleSTD[index]
	
	return (x,y,z),normpdf,AngleSTD

def sigmoidPDF(dist,effLen,distrSlope,sign,highstd=1,loweststd=0):
	### k or [distrSlope] controls the steepness and cellrad*2 controls where the fall will be ###
	return loweststd + (highstd-loweststd)/(1 + np.exp(sign*distrSlope*(dist - effLen)))	

def normalizeVals(x,a=0,b=1):
	return a + ((x-min(x))*(b-a))/(max(x)-min(x))
	
def getRandomOnCellSurface(cellrad,cellcenter,N=1000):
	us = np.random.rand(N)
	vs = np.random.rand(N)
	
	ths = 2*np.pi*us
	phs = np.arccos(2*vs-1)
	
	xs = cellcenter[0] + cellrad*np.sin(phs)*np.cos(ths)
	ys = cellcenter[1] + cellrad*np.sin(phs)*np.sin(ths)
	zs = cellcenter[2] + cellrad*np.cos(phs)

	index = np.random.choice(np.arange(N))
	return xs[index],ys[index],zs[index]
	

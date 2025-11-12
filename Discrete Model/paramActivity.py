import numpy as np

hc = 1.0 # height of the articular cartilage 
hs = 1.0 # height of synovial cavity
hb = 2.0 # bone 

rho = 1.0	# density of synovial fluid
eta = 1.0	# viscosity of synovial fluid
etatc = 0.2	# effective viscosity of synovial fluid insdie the articular cartilage
etatb = 0.2	# effective viscosity of synovial fluid insdie the bone
phic = 0.5 
phib = 0.1

## permeabilities
k = 2.0
k0 = 1.0
k1 = 0.1
kb = 1.0

bet = 0.01		# rate of degradation
alp = 0.01		# rate of deposition
eps = 1.0
D = 0.01
nois = 0.25

Tp = 1.0			     # time period of one gait cycle
T = 300.0*Tp			 # total simulation time
dt = 0.001				 # step size
interval = 10000		 # data saving interval
t = np.linspace(0,T,int(T/dt)+1)

nu = eta/rho
nutc = etatc/rho
nutb = etatb/rho

NbinH = 22
yb,dyb = np.linspace(-hb,0,40,retstep=True)
yc,dyc = np.linspace(0,+hc,NbinH,retstep=True)
ys,dys = np.linspace(hc,hc+hs,20,retstep=True)

#ubx = 0*yb
#ucx = 0*yc
#usx = 0*ys

NbinPhi = 20
ths,dth = np.linspace(-np.pi/2,np.pi/2,NbinPhi,endpoint=False,retstep=True)
Y,Th = np.meshgrid(yc,ths)

" flow params "
overlap = 0.2
us0 = 2.0		## amplitude of extension/flexion of joint
un0 = 0.1		## amplitude of normal loading of joint

" discrete fiber network params "
Wsup,Lsup,Hsup = hc,hc/3,hc		   ## domain size
meanLc 	= Hsup/30				   ## size of the individual fibers
stdLc	= np.round(meanLc*0.2,4)
#discLen = 0.01 #np.round(fdia/0.9,4)
discLen = meanLc
fdia    = round(discLen*0.9,4) 			#discLen/np.sqrt(6) #np.round(cellDia*fib2cell,4) #discLen*0.9 #discLen/np.sqrt(6)	### ref. [Johansson (1993) j. Chem. Phys]	
frad    = 0.5*fdia
#distLc	= np.random.normal(meanLc,stdLc,500) ; 
#distLc  = distLc[distLc > meanLc-3*stdLc] ;	distLc = distLc[distLc < meanLc+3*stdLc]		## limiting the value from [-3*sigma] to [+3*sigma]
#discNp 	= int(meanLc/discLen)
discNp 	= 2
curv    = 0.001
th0		= 0 ; thstd = 5
thetaDist = np.random.normal(th0,thstd*(np.pi/180),1000)

alpha = 5 ; alpstd = alpha*0.2 ; alpdist = np.random.normal(alpha,alpstd,100)			## rate fo fiber deposition in discrete fiber
#tau_beta = np.round(dt*0.9,8)
tau_beta = 0.9998 					## fiber degradation threshold 

Nsample = 10		## number of simulations




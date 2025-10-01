import numpy as np

hc = 1.0 # height of the articular cartilage 
hs = 1.0 # height of synovial cavity
hb = 2.0 # bone 

rho = 1.0
eta = 1.0
etatc = 0.2
etatb = 0.2
phic = 0.5
phib = 0.1

k = 2.0
k0 = 1.0
k1 = 0.1

kb = 1.0

alp = 0.01
eps = 1.0
D = 0.01
bet = 0.01
nois = 0.25

Tp = 1.0 # time period
T = 300.0*Tp
dt = 0.001
interval = 10000
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
us0 = 2.0
un0 = 0.1

" fiber network params "
Wsup,Lsup,Hsup = hc,hc/3,hc
meanLc 	= Hsup/30
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

alpha = 5 ; alpstd = alpha*0.2 ; alpdist = np.random.normal(alpha,alpstd,100)			## rate fo fiber deposition
#tau_beta = np.round(dt*0.9,8)
tau_beta = 0.9998

Nsample = 10




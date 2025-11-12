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
k = 1.0
k0 = 2.0
k1 = 0.1
kb = 1.0

bet = 0.01		# rate of degradation
alp = 0.01		# rate of deposition
eps = 1.0
D = 0.01
nois = 0.25

Tp = 1.0			     # time period of one gait cycle
T = 40.0*Tp				 # total simulation time
dt = 0.001				 # step size
interval = 10000		 # data saving interval

u0 = 1.0

nu = eta/rho
nutc = etatc/rho
nutb = etatb/rho

t = np.linspace(0,T,int(T/dt)+1)
yb,dyb = np.linspace(-hb,0,40,retstep=True)
yc,dyc = np.linspace(0,+hc,20,retstep=True)
ys,dys = np.linspace(hc,hc+hs,20,retstep=True)

ubx = 0*yb
ucx = 0*yc
usx = 0*ys

ths,dth = np.linspace(-np.pi/2,np.pi/2,25,endpoint=False,retstep=True)
Y,Th = np.meshgrid(yc,ths)

delta = 2*np.pi*np.random.rand(1)*0
ux0 = 1.0*u0*np.sin(2*np.pi*t/Tp+delta)#*(t<Tp)
uy0 = 0.1*u0*np.cos(2*np.pi*t/Tp+delta)#*()

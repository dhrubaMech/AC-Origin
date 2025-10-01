import numpy as np
from PIL import Image
import pyvista as pv
import matplotlib.pyplot as plt 

from matplotlib import rc
plt.rcParams["font.family"] = "serif"   
plt.rcParams["mathtext.fontset"] = "stix"
rc('text', usetex=True)
tickFontSize=16
labelFontSize=20
figPanelFontSize=24

def getCMapValues(colname,ncolrs):
	if colname == "hsv":
		colors = plt.cm.hsv(np.linspace(0,1,ncolrs+1))
		#colors = plt.cm.twilight_shifted(np.linspace(0,1,ncolrs+1))
	elif colname == "RdBu":
		colors = plt.cm.RdBu(np.linspace(0,1,ncolrs+1))
	elif colname == "RdBu_r":
		colors = plt.cm.RdBu_r(np.linspace(0,1,ncolrs+1))
	elif colname == "gnuplot":
		colors = plt.cm.gnuplot(np.linspace(0,1,int(ncolrs/2)+1))
		colors = np.concatenate([colors,colors[::-1][1:]])
	elif colname == "custom":
		import matplotlib.colors as mcolors
		#colors = [(0, 'black'), (0.167, 'blueviolet'), (0.333, 'red'), (0.5, 'yellow'), (0.667, 'lime'), (0.8333, 'dodgerblue'), (1.0, 'black')]
		#colors = [(0, 'black'), (0.167, 'blueviolet'), (0.333, 'red'), (0.5, 'yellow'), (0.667, 'darkorange'), (0.8333, 'navy'), (1.0, 'black')]
		#colors = [(0, 'black'), (0.167, 'blueviolet'), (0.333, 'red'), (0.5, 'yellow'), (0.667, 'darkorange'), (0.8333, 'magenta'), (1.0, 'black')]
		#colors = [(0, 'black'), (0.167, 'blueviolet'), (0.333, 'red'), (0.5, 'yellow'), (0.667, 'darkorange'), (0.8333, 'dodgerblue'), (1.0, 'black')]
		colors = [(0, 'black'), (0.167, 'blueviolet'), (0.333, 'red'), (0.5, 'yellow'), (0.667, 'darkorange'), (0.8333, 'cyan'), (1.0, 'black')]
		#colors = [(0, 'black'), (0.167, 'blueviolet'), (0.333, 'red'), (0.5, 'yellow'), (0.667, 'darkorange'), (0.8333, 'teal'), (1.0, 'black')]
		cmap = mcolors.LinearSegmentedColormap.from_list("black_purple_red_yellow_green_blue_black", colors)
		colors = cmap(np.linspace(0,1,ncolrs+1))
		colors = np.concatenate([colors,colors[::-1][1:]])
		#print(colors) ; exit()
	return colors

def plotfibers(xf0,zf0,ph0fibs,Nfibs,colors,colIdx,ni,Wsub,Hsub,phFibdist,kxx,kyy,ax):
	ax[0].clear()
	image = Image.open("fib.png")
	ax[0].imshow(image,zorder=0)
	ax[0].set_xlim([125,673])
	ax[0].set_xticks([125,400,673],["$0$","$0.5$","$1.0$"])
	ax[0].set_ylim([675,125])
	ax[0].set_yticks([675,400,125],["$0$","$0.5$","$1.0$"])
	
	ax[0].set_title("$t_i=$"+f" {ni} | "+"$N_{fibs}$"+f" : {Nfibs}",fontsize=tickFontSize)
	ax[0].set_ylabel("$z_s$",fontsize=labelFontSize) ; ax[0].set_xlabel("$x_s$",fontsize=labelFontSize)
	
	ax[1].clear()
	ax[1].contourf(Th,Y,phFibdist,200)
	ax[1].set_xlabel("$\phi$",fontsize=labelFontSize)
	ax[1].set_ylabel("$z_s$",fontsize=labelFontSize)
	ax[1].set_xticks([-np.pi/2,ths[int(NbinPhi*0.5)-1]+dth*0.5,ths[-1]],["$-\\pi/2$","$0$","$\\pi/2$"])
	
	ax[2].clear()
	ndist = 5 ; 
	#cols = ["r" , "g" , "b" , "orange" , "magenta" , "cyan"]
	cols = plt.cm.jet(np.linspace(0,0.9,ndist))
	idz = np.linspace(1,NbinH,ndist,dtype=int)-1
	for ii in range(ndist):
		ax[2].plot(ths,phFibdist[:,idz[ii]]/np.max(phFibdist),c=cols[ii])
	ax[2].set_xlabel("$\phi$",fontsize=labelFontSize)
	ax[2].set_ylabel("$p_\\phi$",fontsize=labelFontSize)
	ax[2].set_xticks([-np.pi/2,0,+np.pi/2],["$-\\pi/2$","$0$","$\\pi/2$"])
	ax[2].set_ylim([0,1.0])
	
	ax[3].clear()
	ax[3].plot([0,0],[0,1],c="gray",ls="--",lw=0.5)
	#ax[3].plot(1/kxx - 1/kyy,yc,'k')
	ax[3].plot((kxx-kyy)/(kxx+kyy),yc,'r')
	ax[3].plot([-0.045,-0.045],[0,1],c="gray",ls="--",lw=0.5)
	ax[3].plot([+0.045,+0.045],[0,1],c="gray",ls="--",lw=0.5)
	ax[3].set_ylabel("$z_s$",fontsize=labelFontSize)
	ax[3].set_xlabel("$\\kappa$",fontsize=labelFontSize)
	#ax[3].set_xlim([-0.5,0.5])
	ax[3].set_xlim([-0.11,0.06])
	ax[3].set_ylim([0,1])
	
	for i in range(4):
		ax[i].tick_params(labelsize=tickFontSize)
	
	plt.tight_layout()
	plt.draw()
	plt.pause(1e-10)
	#print(input())

def pltDist(ucy,ucx,p,kxx,kyy):
	phi = np.arctan2(ucy,ucx)
	phi[phi < -np.pi/2] += np.pi
	phi[phi >  np.pi/2] -= np.pi

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

	plt.pause(0.001)
	ax[0,0].cla()
	ax[0,1].cla()
	ax[1,0].cla()
	ax[1,1].cla()

def renderfiber(pl,fibdata,xf0,yf0,zf0,Nfibs,ph0fibs,th0fibs,colors,colIdx,meanLc,name):
	
	x2 = xf0 + meanLc*np.cos(ph0fibs)
	y2 = yf0 + meanLc*np.sin(th0fibs)
	z2 = zf0 + meanLc*np.sin(ph0fibs)
	xf = np.column_stack((xf0,x2))
	yf = np.column_stack((yf0,y2))
	zf = np.column_stack((zf0,z2))
	xf += np.repeat(xf0 - np.mean(xf,axis=1),2).reshape(Nfibs,2)
	yf += np.repeat(yf0 - np.mean(yf,axis=1),2).reshape(Nfibs,2)
	zf += np.repeat(zf0 - np.mean(zf,axis=1),2).reshape(Nfibs,2)
	
	fibs = np.zeros(3*Nfibs,dtype=int)
	colfibs = []
	for fi in range(Nfibs):
		fibs[fi*3] = 2 
		fibs[fi*3+1] = fi*2 
		fibs[fi*3+2] = fi*2+1 
		colfibs.append(colors[int(colIdx[fi])][:3])
	
	fibdata.points = np.column_stack((xf.flatten(),yf.flatten(),zf.flatten()))
	fibdata.lines = fibs.flatten()
	fibdata.cell_data['Orientation'] = np.array(colfibs)
	sargs = dict(height=0.25, vertical=True, position_x=0.9, position_y=0.40)
	pl.add_mesh(fibdata,color='white',render_lines_as_tubes=True,line_width=5.0,smooth_shading=True,ambient=0.5,diffuse=1.0,specular=0.0,scalars='Orientation',scalar_bar_args=sargs,rgb=True)
	#pl.screenshot(f'{name}.png',window_size=[800,800])
	pl.view_xz()
	pl.enable_parallel_projection()
	pl.screenshot(f'fib.png',window_size=[800,800]) 

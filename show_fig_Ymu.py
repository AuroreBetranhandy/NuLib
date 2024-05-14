import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib.colors import LogNorm, SymLogNorm

import h5py

init_data= h5py.File('SRO.h5','r')


EOS_mu_mu = np.asarray(init_data["mu_mu"])
EOS_mu_e = np.asarray(init_data["mu_e"])
Ye = np.asarray(init_data["ye"])
EOS_temp = np.asarray(init_data["logtemp"])
EOS_dens = np.asarray(init_data["logrho"])
EOS_temppoint = list(init_data["pointstemp"])
EOS_denspoint = list(init_data["pointsrho"])

init_data.close()



diff_mu = np.load("diff_muon_full.npy")
result_mu = np.load("result_muon_full.npy")

diff_e = np.load("diff_elec_full.npy")
result_e = np.load("result_elec_full.npy")



	# ~ fig,ax = plt.subplots(1,3,figsize=(15.,10.))


print(Ye[30])


Ye_temp= Ye[10]
T_temp=EOS_temp[20]
rho_temp=EOS_dens[100]


print(Ye_temp,result_mu[10,20,100])
print(np.where(abs(Ye -result_mu[10,20,100])< 1e-2)[0][0])
loc=np.where(abs(Ye -result_mu[10,20,100])< 1e-2)[0][0]

print(EOS_mu_mu[loc,20,100],EOS_mu_mu[10,20,100])
# ~ C = np.log10(EOS_mu_e[30,:,:])
# ~ C[C<-3]=np.nan
# ~ result_mu[result_mu>1.]=np.nan
# ~ x,y,z=np.meshgrid(Ye,EOS_temp,EOS_dens)
# ~ fig = plt.figure()
# ~ ax = fig.add_subplot(111, projection='3d')
# ~ im = ax.scatter(x,y,z,c=result_mu.T)
# ~ fig.colorbar(im)
# ~ plt.show()
for i in range(0,51,10):
	fig,ax = plt.subplots(1,3,figsize=(15.,10.))
	ax1=ax[0]
	ax2=ax[1]
	ax3=ax[2]
	im1 = ax1.pcolor(EOS_dens,EOS_temp,result_mu[i],cmap='jet',vmin=0,vmax=0.51)
	divider = make_axes_locatable(ax1)
	cax = divider.append_axes('left', size='5%', pad=1.)
	fig.colorbar(im1,cax=cax, orientation='vertical')
	
	im2 = ax2.pcolor(EOS_dens,EOS_temp,(diff_mu[i]),cmap='seismic',vmin=-1,vmax=1)
	divider = make_axes_locatable(ax2)
	cax = divider.append_axes('left', size='5%', pad=1.)
	fig.colorbar(im2,cax=cax, orientation='vertical')
	
	im3 = ax3.pcolor(EOS_dens,EOS_temp,(result_e[i]-result_mu[i]),cmap='jet',vmin=0,vmax=0.51)
	divider = make_axes_locatable(ax3)
	cax = divider.append_axes('left', size='5%', pad=1.)
	fig.colorbar(im3,cax=cax, orientation='vertical')

	fig.suptitle(str(Ye[i]))
	plt.show()


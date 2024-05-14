import numpy as np
import matplotlib.pylab as plt

energies=np.loadtxt("muons.txt", unpack=True,usecols=(0))
scat_e= np.loadtxt("muons.txt", unpack=True,usecols=(1,2,3,4,5,6))
scat_mu= np.loadtxt("muons.txt", unpack=True,usecols=(7,8,9,10,11,12))


fig,ax = plt.subplots(3,2,figsize=(20,10))

species=[ "nue","anue","numu","anumu","nutau","anutau"]

for i in range(6) : 
	ax[abs((i-3)%3),int(i/3)].plot(energies,scat_e[i,:],label="e")
	ax[abs((i-3)%3),int(i/3)].plot(energies,scat_mu[i,:],label='mu')
	ax[abs((i-3)%3),int(i/3)].set_title(species[i])
	ax[abs((i-3)%3),int(i/3)].set_ylim([0.1,100])
	ax[abs((i-3)%3),int(i/3)].legend()
plt.show()




energies_em=np.loadtxt("muons_emi.txt", unpack=True,usecols=(0))
emi= np.loadtxt("muons_emi.txt", unpack=True,usecols=(1,2,3))
abso= np.loadtxt("muons_emi.txt", unpack=True,usecols=(4,5,6))


fig,ax = plt.subplots(1,2,figsize=(20,10))

ax[0].semilogy(energies_em,emi[0,:],'.-b')
ax[0].semilogy(energies_em,emi[1,:],'r')
ax[0].semilogy(energies_em,abso[1,:],':r')

ax[1].semilogy(energies_em,emi[2,:])
ax[1].semilogy(energies_em,abso[2,:],':')

ax[0].set_ylim([1e-4,1.])
ax[1].set_ylim([1e-5,1e-1])

plt.show()



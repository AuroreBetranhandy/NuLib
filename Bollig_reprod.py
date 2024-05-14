import numpy as np

import h5py 
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

import numpy as np
import matplotlib.pylab as plt
import sys
from scipy import interpolate



sys.path.insert(1, '/home/aube8745/PhD/NuLib-muons/eos')
sys.path.insert(2, '/home/aube8745/PhD/NuLib-muons/helmotz')
import NuclearEos as ne
import helmholtz
table="SFHo_low.h5"
#table="SRO.h5"
neos = ne.NuclearEOS(table)

# ~ sys.exit()

#### mass excesses in MeV
u = 931.49410242

m_n_excess = 0.008664*u
m_p_excess = 0.007276*u
m_alpha_excess =  0.001506179127*u
m_iron_excess = -60.603



### masses in MeV 
mev_to_gram = 1.782661758e-27

m_n = 939.565346 
m_p = 938.272013
m_alpha = 3727.379240
m_iron = 60 *u #55.845


m_mu=105.66
m_e=0.511
hbarc_mevcm = 1.97326966e-11# in MeV*cm



kelvin_to_mev = 8.6173423e-11
G = 6.674e-8 * 2e33**2
number_of_zones = 796
c_light = 3.e10
mev_to_erg = 1.60217733e-6


var = ne.EOSVariable()
var.xtemp = 0.01
var.xrho = 1e10
var.xye=0.1
#var.xye = Hann_values[200,4,0]
var = neos.nuc_eos_full(var,mode=ne.EOSMODE_RHOT)

print(var.xmu_e)


# ~ xmass=np.zeros([4]) ; aion=np.zeros([4]) ; zion=np.zeros([4]) 

# ~ xmass[0] = var.xxn ; aion[0]  = 1.0 ; zion[0]  = 0.0
# ~ xmass[1] = var.xxp ; aion[1]  = 1.0  ; zion[1]  = 1.0
# ~ xmass[2] = var.xxa ; aion[2]  = 4.0  ; zion[2]  = 2.0
# ~ xmass[3] = var.xxh ; aion[3]  = 30.583 ; zion[3]  = 14.5
# ~ #~ var.xrho = [Hann_values[200,1,5],Hann_values[200,1,66]]




x = np.linspace(0,500,500)

T=np.logspace(-3,2,100)
rho_ye=np.logspace(-2,15,200)





def int_energ(T, m_l, mu,x):
	e_int = 8. * np.pi * T**4 /(2. * np.pi * hbarc_mevcm)**3 
	integral = np.sum((x+ m_l/T) * np.sqrt(x*(x + 2 * m_l /T) / (np.exp( x - (mu - m_l)/T) + 1)))
	return e_int*integral

matrix = np.zeros((100,200))
for i in (range(100)):
	for j in range(200):
		var.xtemp = T[i]
		var.xrho = rho_ye[j]
		var.xye=0.5
		var = neos.nuc_eos_full(var,mode=ne.EOSMODE_RHOT)

		matrix[i,j] = var.xmu_e #int_energ(var.xtemp,m_mu,var.xmu_e,x)


fig, ax = plt.subplots()

ax.set_yscale('log')
ax.set_xscale('log')
ax.pcolor(rho_ye,T, np.log10(matrix))



#plt.imshow(matrix, extent=(T[0],rho_ye[0],T[-1],rho_ye[-1]),aspect='auto')
plt.show()

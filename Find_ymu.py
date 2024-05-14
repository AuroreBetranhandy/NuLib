import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from scipy import integrate

import h5py

init_data= h5py.File('SRO.h5','r')


EOS_mu_mu = np.asarray(init_data["mu_mu"])
EOS_mu_e = np.asarray(init_data["mu_e"])
Ye = list(init_data["ye"])
EOS_temp = list(init_data["logtemp"])
EOS_dens = list(init_data["logrho"])
EOS_temppoint = list(init_data["pointstemp"])
EOS_denspoint = list(init_data["pointsrho"])

init_data.close()


print(len(Ye))
print((EOS_temp[-1]))


clight=3.e10 #cm/s
hbarc_mevcm = 1.97326966e-11 # in MeV*cm
m_mu = 105.66  #MeV
# ~ m_mu = 0.511 #105.66  #MeV
avo = 6.0221367e23
m_amu =  931.494061

mev_to_gram = 1.782661758e-27 #one MeV is # grams
def integrand(beta,xi):
	nyumph=0.
	integrand_space=np.logspace(-4,5,num=1000)
	integrand_dx=integrand_space[1:]-integrand_space[:-1]
	for i in range(1,len(integrand_space)-1):
		x = integrand_space[i]
		dx= integrand_dx[i]
		nyumph = nyumph+ ((x+beta)*np.sqrt(x*(x+2.*beta))/(np.exp(x-xi)+1.))*dx
		
	return nyumph
	





def find_ymu(rho,T,mu_mu):
	 constant= 8.*np.pi*T**3 /(hbarc_mevcm *2.*np.pi )**3	
	 beta=m_mu/T
	 eta=(mu_mu)/T
	 
	 xi_min= eta - beta
	 xi_plus = -xi_min - 2*(m_mu/T)
	 inte_min= integrand(beta,xi_min)
	 inte_plus= integrand(beta,xi_plus)
	 #print(inte_min,inte_plus,constant,eta,beta)
	 return  ((inte_min-inte_plus)*constant)  /(rho/(m_amu*mev_to_gram))

Temp = 10**(EOS_temp[132])
rho = 10**(EOS_dens[200])
mu_mu = EOS_mu_mu[30,132,200]

num_rho= len(EOS_dens)
num_temp= len(EOS_temp)
num_mu= len(Ye)
result= np.zeros((num_mu,num_temp,num_rho))
diff= np.zeros((num_mu,num_temp,num_rho))
Y=np.zeros((num_mu,num_temp,num_rho))

step=1
for  k in range(0,num_mu,step):
	print(k/num_mu)
	for i in range(0,num_rho,step): 
		for j in range(0,num_temp,step): 
			Temp = 10**(EOS_temp[j])
			rho = 10**(EOS_dens[i])
			mu_mu = EOS_mu_mu[k,j,i]
			result[k,j,i]=find_ymu(rho,Temp,mu_mu) 
			diff[k,j,i]=abs(result[k,j,i]-Ye[k])

np.save("diff_muon_full",diff)
np.save("result_muon_full",result)



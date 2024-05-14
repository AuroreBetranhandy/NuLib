import numpy as np
import matplotlib.pylab as plt
import sys
from scipy import interpolate
import pickle 
import multiprocessing as mp
import h5py


init_data= h5py.File('SRO.h5','r')

#print(list(init_data["mu_e"]))

Y_mu = np.load('result_muon_full.npy')
Y_e = np.load('result_elec_full.npy')

init_data.close()


EOS = h5py.File('SRO.h5','a')



EOS.create_dataset('Y_e_new', data=Y_e-Y_mu)

EOS.close()

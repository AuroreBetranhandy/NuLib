import csv
import numpy as np
import h5py 
# n (fm-3) 	ε (MeV/fm3) 	P (MeV/fm3) 	x p 	xe

m_amu =  931.494061 #MeV
Avog_num=6.022e23

test_matrix=np.loadtxt('EOS_muons.csv', unpack=True)

density= test_matrix[0]*m_amu/Avog_num*(1e12)**3


hf = h5py.File('EOS.h5', 'w')
dset = hf.create_dataset("rho", data=density)
dset = hf.create_dataset("energ", data=test_matrix[1])
dset = hf.create_dataset("P", data=test_matrix[2])
dset = hf.create_dataset("x_p", data=test_matrix[3])
dset = hf.create_dataset("x_e", data=test_matrix[4])

import h5py
import numpy
from numpy import *
from scipy import *
from sys import *

f=h5py.File(argv[1])

ds_freq = f["frequencies"] 
ds_sqw = f["sqw"]

freq = ds_freq[:]
sqw = ds_sqw[:]

for i in range(len(sqw)):
 print freq[i],sqw[i,0], sqw[i,1]

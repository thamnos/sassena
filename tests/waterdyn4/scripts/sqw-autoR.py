import h5py
import numpy
from numpy import *
from scipy import *
from sys import *

def gaussian(x, sigma = 1.0, mu = 0.0, normalize = False):
    """Returns a gaussian.
    """
            
    gaussian = exp(-0.5*( (x - mu) / sigma)**2)
    
    if normalize:
    
        gaussian *= 1.0/(sigma*sqrt(2.0*pi))

    return gaussian
    
def gaussianWindow(inputSeries, gaussian):
    """Returns a gaussian smoothed signal.
    
    @param inputSeries: the signal to smooth.
    @type inputSeries: NumPy array   

    @param alpha: a float specifying the width of the smoothing gaussian.
    @type alpha: float

    @return: the smoothed signal.
    @rtype: NumPy array
    """

    # outputSeries is an array of length 2*len(series)-1 directionized to zero 
    # outputSeries is the smoothed version of inputSeries obtained by applying
    # a gaussian kernel to intputSeries
    outputSeries = zeros((2*len(inputSeries) - 2,), dtype = complex)        
        
    # exp(...) is the gaussian window used to smooth the spectrum series 
    res = inputSeries*gaussian

    # The second half of outputSeries is filled using periodic conditions
    outputSeries[:len(inputSeries)] = res
    outputSeries[len(inputSeries):] = res[-2:0:-1]

    return outputSeries

#print "Got a total of ", len(argv), " arguments (including script file)"
if (len(argv)<4):
    print "need 3/5 arguments: hdf5 signal file, qindex, hdf5 output file, [ fwhm_e , dt ]"
    exit(1)

f = h5py.File(argv[1])
if (not "fqt" in f):
    print "hdf5 signal file doesn't contain the dataset fqt"
    exit(1)
    
qindex = int(argv[2])

ds = f["fqt"]

fqt = zeros(ds.shape[1],dtype=complex)
fqt.real = ds[qindex,:,0]
fqt.imag = ds[qindex,:,1]

if (len(argv)<5):
    fwhm_t = float(int(0.1*len(fqt)))
    sigma_t = 2*fwhm_t / ( 2*sqrt(2*log(2)) )
    dt = 1
else:
    fhmw_e = float(argv[4]) # FHMW in meV
    dt = float(argv[5]) # in ps
    sigma_e = fhmw_e / 2.354820045
    sigma_t = 1.0/(1.5192669*sigma_e)

time = arange(len(fqt))*dt
frequencies = arange(len(time))/(2.0*len(time)*dt)
r_t = gaussian(time,sigma_t,0)


fqtgaussian = gaussianWindow(fqt,r_t)
sqw = fft(fqtgaussian)

fout = h5py.File(argv[3],"w")
dsout = fout.create_dataset("sqw", ((len(sqw),2)),"=f8")
dsout2 = fout.create_dataset("frequencies",((len(sqw),)))

dsout[:,0]=sqw[:(len(sqw))].real
dsout[:,1]=sqw[:(len(sqw))].imag
dsout2[:(len(sqw)/2)]=frequencies[:(len(sqw)/2)]
dsout2[(len(sqw)/2):]=-1*frequencies[(len(sqw)/2):0:-1]

fout.close()

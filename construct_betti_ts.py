import chomp_betti as CB
import rbc_basic as R
import os
import matplotlib.pyplot as plt
from pylab import *

def make_timeseries( files, length=-1 ):

    frames= os.listdir( files )
    bettis = []
    
    frames.sort( key=R.natural_key )    
    for frame in frames[:length]:
        bettis.append( loadtxt( files+frame ) )
    arr = asarray( bettis, dtype=int )
    return arr

def plot_timeseries( arr, dim=0, cell_name="", fft_pass=0.05 ):
    """
    arr.shape = ( frame #, dimension, generators )
    """
    bettis = arr[:,1]

    fig = figure()
    ax = fig.gca()
    ax.set_title( r"$H_{"+str(dim)+"}$ time series for cell "+ cell_name +\
                 " with lowpass of "+str(fft_pass) )
    ax.plot( bettis )
    #fig.savefig( savedir+'betti_deriv_dim'+str(dim) )
    

if __name__ == "__main__":

    #files = '/data/rbc_deriv_files/RBC_Deriv_110125_r01/'
    files = '/Users/jberwald/Dropbox/Projects/rbc/RBC_Wyss/rbc_shared/betti_ts/'
    cell = 'RBC_Deriv_110125_r005/'
    cell_name = 'Deriv_110125'
    #subdir = 'betti/'
    savedir = files

    # script to concatentate a bunch of separate chomp output files
    # into one array
    #bt = make_timeseries( length=-1 )
    #bt0 = bt[:,0,:]
    #bt1 = bt[:,1,:]
    #savetxt( files+'betti_ts_dim0', bt0 )
    #savetxt( files+'betti_ts_dim1', bt1 )

    # script to plot the betti time series
    dim = 1
    bt = loadtxt( files+cell+'betti_ts_dim'+str(dim) )
    plot_timeseries( bt, dim=dim, cell_name=cell_name )

    show()

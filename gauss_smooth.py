from scipy import *
from scipy import signal
import cPickle as pkl
import numpy as np

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return improc

if __name__ == "__main__":

    import sys

    new_fdir = '/sciclone/home04/jberwald/data10/wyss/data/Cells_Jesse/New/frames/'
    old_fdir = '/sciclone/home04/jberwald/data10/wyss/data/Cells_Jesse/Old/frames/'
    new_cells = [ 'new_110125/', 'new_130125/', 'new_140125/', 'new_40125/', 'new_50125/' ]
    old_cells = [ 'old_100125/', 'old_120125/', 'old_50125/', 'old_90125/' ]

    #  modes = [ ('04','05'), ('05','06'), ('06','07'),('07','08'),('08','09')]

    arg = int( sys.argv[1] )
    #modes_arg= int( sys.argv[2] )
    #modes
    rlow = '07' #modes[0]
    rhigh = '08' #modes[1]
    path = new_fdir + new_cells[ arg ] + 'fft_frames/bandpass/'
    fname = 'fft_r'+rlow+'_r'+rhigh+'.pkl'

    print "Loading data from", path+fname
    with open( path+fname ) as fh:
        data = pkl.load( fh )
    # smooth the images
    print "smoothing images..."
    frames = [ frame['ifft_nonzero'] for frame in data ]
    smoothed = [ blur_image( f, 3 ) for f in frames[:1000] ]
    masked = [ np.ma.masked_less( f, f.mean() ).mask for f in smoothed ]
    savename = path + 'smooth_r'+rlow+rhigh+'/smooth_mask_'
    print "saving smoothed images..."
    print ""
    for i,ma in enumerate( masked ):
        np.save( savename+str(i), ma )

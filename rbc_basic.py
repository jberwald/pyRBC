"""
Module containing some basic extraction functions to obtains cell frames in a list (load_rbc()).

Useage: ('In [m]' refers to input number m in the ipython interpreter.)

In [1]: fname = '/data/jberwald/rbc/New/new_140125-concatenated-ASCII'

( Replace fname with the path to a cell on the whatever system you're
on. Eg., on simplex

fname = 'data/jberwald/wyss/data/Cells_Jesse/New/cells/new_140125-concatenated-ASCII'

should work. )

In [2]: import rbc_basic as R

In [3]: f = R.load_rbc( fname, 4997, 203,198 )

In [4]: R.plot_frame( f[0] )

(You should see the image 'rbc_basic.png' included in the same
directory as this module. Note the halo--the boundary has not been
removed.)

"""

import numpy
#from pylab import *
import re
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import rbc_processing as rp
import sys

slash = '/'

def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def load_rbc( fname, skiprows, nx, ny ):
    """
    Returns a block of frames from <fname> cell data. Reshapes array
    according to boundary data determined from nx, ny (determined from
    corresponding boundary file).
    """
    C = numpy.loadtxt( fname, skiprows=skiprows )    
    cell_frames = [ C[i].reshape(( nx,ny )) for i in range( 5000-skiprows ) ]
    return cell_frames

def extract_frames( fname, bnd ):
    """
    Each line of the cell in <fname> contains a raveled matrix. Read
    each line of <fname>, reshape it based on array stored in <bnd>
    file. Then save each frame to file.
    """
    if os.uname()[0] == 'Linux':
        savefunc = numpy.save
    elif os.uname()[0] == 'Darwin':
        # why the !@$# doesn't Mac save readable .npy files??
        savefunc = numpy.savetxt
    fh = open( fname, 'r' )
    bnd_arr = numpy.loadtxt( bnd )
    nx, ny = bnd_arr.shape
    bnd_arr = bnd_arr.ravel()
    # save files in a special folder
    part = fname.rpartition( '/' )
    cell_name = part[-1]
    # find the subfolder name from the beginning of the cell
    # name. subfolder must have been created already.
    subfolder = cell_name.partition( '-' )[0]
    savedir = part[0] + '/' + 'frames/' + subfolder +'/'

    make_dir( savedir )
    # print 'savedir', savedir
     
    # loop over lines in <fname> and save each as nx x ny array
    #k = 0. 
    fromstring = numpy.fromstring
    for k, line in enumerate( fh.readlines() ):
        arr = fromstring( line, sep='\t' )
        # remove boundary 
        arr = bnd_arr * arr
        arr.resize( (nx,ny) )
        #numpy.save( savedir + cell_name + '_' + str(k), arr )
        savefunc( savedir + cell_name + '_' + str(k), arr )

def plot_frame( frame ):
    """
    Use pylab.imshow() to plot a frame. Note: if the bounday has not
    been extracted, then there will be a halo. (See # remove bounday
    in extract_frames() above.)
    """
    fig = plt.figure()
    ax = fig.gca()
    ax.set_title( "RBC frame" )
    ax.imshow( frame )
    fig.show()

def plot_frame_mask_zero( frame, nx=203, ny=198 ):
    """
    Use pylab.imshow() to plot a frame. Mask the elements of the image
    matrix that are zero.
    """
    frame.resize((nx,ny))

    cdict = {'red': ((0., 1, 1),
                     (0.05, 1, 1),
                     (0.11, 0, 0),
                     (0.66, 1, 1),
                     (0.89, 1, 1),
                     (1, 0.5, 0.5)),
            'green': ((0., 1, 1),
                      (0.05, 1, 1),
                      (0.11, 0, 0),
                      (0.375, 1, 1),
                      (0.64, 1, 1),
                      (0.91, 0, 0),
                      (1, 0, 0)),
            'blue': ((0., 1, 1),
                     (0.05, 1, 1),
                     (0.11, 1, 1),
                     (0.34, 1, 1),
                     (0.65, 0, 0),
                     (1, 0, 0))}

    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
    fig = plt.figure( frameon=False, dpi=160 )
    ax = fig.gca()
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    im = ax.imshow( frame, cmap=my_cmap )
    cbar = plt.colorbar( im, shrink=0.6 )
    fig.show()
    return fig
    
def plot_sublevel_set( fname, bndfile=None, persfile=None,
                       height=None, mask=False, nx=203, ny=198 ):
    """
    Plot sublevel set for a cell <fname> (full path to cell).
    """
    data = numpy.loadtxt( fname )
    print data.shape
    if bndfile:
        bnd = numpy.loadtxt( bndfile )
        print bnd.shape
        data = bnd.ravel() * data
        data.resize( (nx,ny) )
    if persfile:
        heightList = [x[1:-1] for x in re.findall('\[[^\]]*\]',s)]
    elif height:
        heightList = [height]
    else:
        print "Must provide either a path to a Perseus output file or a single height for the sublevel set"
        sys.exit( 1 )

    # make an RBG array to hold (x,y,z) data at each pixel
    G = numpy.zeros((nx,ny,3), dtype=int)
    
    outdir = slash.join( fname.split( '/' )[:-1] ) + '/'
    for ind, h in enumerate(heightList):
        #temp = data.copy()
        G[ numpy.where( data > int(h) ) ] = [1,1,1]
        if mask:
            G[ numpy.where( data <= int(h) ) ] = [0,0,160]
            G[ numpy.where( data == 0 ) ] = [1,1,1]
        outName = fname.split('/')[-1][:-4] + '_' + str( h )
        output = outdir + outName
        # now plot stuff
        fig = plt.figure( frameon=False )
        ax = fig.gca()
        #ax.set_title('RBC ' + output)
        ax.imshow( G )
        ax.set_xticks( [] )
        ax.set_yticks( [] )
        fig.show()
        #plt.colorbar()
        print "saving images to", output
        fig.savefig( output + '.png', dpi=160 )
        fig.savefig( output + '.pdf', dpi=160 )
    return G


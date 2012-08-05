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
from pylab import *
import re
import sys


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

def cell2npy ( cell_file, bnd_file, shape=(182,196) ):
    """
    Convert a cell file with one frame per line to a numpy file
    containing a 3d array of matrices.
    """
    # unroll the shape tuple
    rows, cols = shape
    bnd = numpy.loadtxt ( bnd_file, dtype=int )
    
    if len( bnd.shape ) > 1:
        bnd = bnd.ravel()

    print "bnd.shape", bnd.shape[0]

    if bnd.shape[0] != rows*cols:
        print "WRONG SHAPE!!"
        sys.exit(1)
    
    # deal with windows' new line issues
    frames = []
    print "Converting txt file to numpy arrays..."
    with open ( cell_file, 'rU' ) as fh:
        for line in fh.readlines():
            tmp = bnd * numpy.fromstring( line, sep='\t', dtype=int )
            #tmp = numpy.fromstring( line, sep='\t', dtype=int )
            frames.append ( tmp.reshape ( shape ) )
    print "Stacking frames..."
    cell = numpy.dstack ( frames )
    del frames # clear some memory 

               #return frames[0], bnd

    print "Saving file..."
    numpy.save ( cell_file+'.npy', cell )
    print "Saved file:", cell_file+'.npy'

def plot_frame( frame ):
    """
    Use pylab.imshow() to plot a frame. Note: if the bounday has not
    been extracted, then there will be a halo. (See # remove bounday
    in extract_frames() above.)
    """
    fig = figure()
    ax = fig.gca()
    ax.set_title( "RBC frame" )
    ax.imshow( frame )
    fig.show()
    

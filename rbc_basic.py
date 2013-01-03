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
import re
#import matplotlib
#matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
#import rbc_processing as rp
import sys, os

slash = '/'

def cell2npy( cell, bnd, rows=(203,198) ):
    """
    Convert a concatenated text file to a NPY file.

    fname -- full path to concatenated file

    bname -- full path to boudnary file

    row -- dimensions (tuple)
    """
    bnd_arr = numpy.loadtxt( bnd )
    frames = []
    print "opening ", cell
    with open( cell, 'rU' ) as fh:
        fromstring = numpy.fromstring
        print "converting lines..."
        for line in fh.readlines():
            arr = fromstring( line, sep='\t' )
            # remove boundary 
            arr = bnd_arr.ravel() * arr
            arr.resize( rows )    
            frames.append( arr )
    savefile = cell + '.npy'
    print "saving file to", savefile
    stack = numpy.dstack( frames )
    numpy.save( savefile, stack )
    return stack
    

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

def extract_frames( fname, bnd, to_png=True, save_prefix='/data/jberwald/old_png/old12' ):
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

    # make_dir( savedir )
    # print 'savedir', savedir
     
    # loop over lines in <fname> and save each as nx x ny array
    #k = 0.
    # ma
    if to_png:
        fig = plt.figure()
    fromstring = numpy.fromstring
    for k, line in enumerate( fh.readlines() ):
        arr = fromstring( line, sep='\t' )
        # remove boundary 
        #arr = bnd_arr * arr
        arr.resize( (nx,ny) )
        if to_png:
            ax = fig.gca()
            im = ax.imshow( arr )
            fig.savefig( save_prefix + 'noBND_' + str( k ) + '.png' )
            fig.clf()
        if k==10:
            return
        #numpy.save( savedir + cell_name + '_' + str(k), arr )
        #savefunc( savedir + cell_name + '_' + str(k), arr )


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
    # This allows frame to be any data, such as a symmetric 2D
    # Gaussian. 
    try:
        frame.resize((nx,ny))
    except:
        pass

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
    # for transparency
    fig.patch.set_alpha( 0.0 )
    ax = fig.gca()
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    ax.set_axis_off() # turn off the axes 
    im = ax.imshow( frame, cmap=my_cmap )
    cbar = plt.colorbar( im, shrink=0.6 )
    fig.show()
    return fig
    
def plot_sublevel_set( frame, height, bndfile=None, persfile=None,
                       nx=203, ny=198, save=False,
                       transparent=True, thresh=None ):
    """
    Plot sublevel set for an array representing an intensity function.

    frame -- path to array on disk or Numpy array.

    height -- function height at which to take sublevel set.

    Returns figure object
    """
    h = height
    # text file
    try:
        try:
            data = numpy.loadtxt( frame )
            # numpy file
        except ValueError:
            data = numpy.load( frame )
    except:
        # This is for a general array that is already rectangular
        # and needs no help.
        data = frame
        nx,ny = data.shape

    print "data", data.shape
    # if we read in a data file we might need a boundary file to go with it.
    if bndfile:
        bnd = numpy.loadtxt( bndfile )
        print "boundary", bnd.shape
        try:
            data = bnd.ravel() * data
        except ValueError:
            data = bnd * data
        data.resize( (nx,ny) )

    # make an array to hold (R,G,B,A) data at each pixel
    G = numpy.zeros((nx,ny,4), dtype=int)

    # if not thresh:
    #     thresh = 0.9
    
    #temp = data.copy()
    G[ numpy.where( data > int(h) ) ] = [1,1,1,0]
    # the sublevel set
    G[ numpy.where( data <= int(h) ) ] = [0,0,160,1]
    # everything outside 
    G[ numpy.where( data == 0.0 ) ] = [1,1,1,0]
    # outside majority of Gaussian peak
    if thresh:
        G[ numpy.where( data <= thresh ) ] = [1,1,1,0]
   
    # now plot stuff
    fig = plt.figure( frameon=False, dpi=160 )

    # make things transparent?
    if transparent:
        fig.patch.set_alpha( 0.0 )

    ax = fig.gca()
    ax.set_frame_on( False )
    #ax.set_title('RBC ' + output)
    # PLOT THE MATRIX OF VALUES
    print "plotting sublevel set..."
    ax.imshow( G )
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    fig.show()
    #plt.colorbar()
    if save:
        if type( frame ) == str:
            # make output directory
            outdir = slash.join( frame.split( '/' )[:-1] ) + '/'
            outName = fname.split('/')[-1][:-4] + '_' + str( h )
            output = outdir + outName
    
        print "saving images to", output
        fig.savefig( output + '.png', dpi=160 )
        fig.savefig( output + '.pdf', dpi=160 )
    return G, fig

def sublevel_sequence():
    heights = [ 900,1050,1110 ]
    frame = 'new11_frame2000.txt'
    bnd = 'boundary_Nov_new110125'
    figs = []
    for h in heights:
        figs.append( plot_sublevel_set( frame, h, bndfile=bnd, mask=True ) )
    return figs


if __name__ == "__main__":

    pass

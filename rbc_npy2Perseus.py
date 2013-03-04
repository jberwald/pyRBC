"""
    
Module to convert RBC frames (npy files) to input required by Perseus software (for computing topological persistence)
created by kel 5/23/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
    
    
    NOTE: There's been found to be discrepancies between sparse and dense cubical
    formats
    - Sparse returns persistence intervals more in line with intuition
"""

import numpy
import matplotlib.pylab as plt
import re
import os
import pickle as pkl
from scipy.spatial import distance


def write_file ( fname, output ):
    """ 
        ----DENSE FORMAT---
        - Write contents of single frame to output file
            e.g. new cell name on simplex:
            '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/new_110125-concatenated-ASCII_324.npy'
            e.g. old cell name on simplex:
            '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_100125/old_100125-concatenated-ASCII_1525.npy'
        Perseus input requires following format:
        #Dimension of cubical grid
        #Number of cubes in first dimension
        #Number of cubes in second dimension
        #Number of cubes in   ...     dimension
        #Birth time of cube anchored at (0,0,...,0)
        #Birth time of cube anchored at (0,0,...,1)
        #Note: Output is file name
    """
    #Load npy file
    data = numpy.load(fname)
    size = data.shape#(203,198)
    dim = len ( size )
    #Open file
    text = open (output + ".txt", "w")
    #Write number of dimensions
    text . write (str(dim) + "\n")
    #Write resolutions
    for i in xrange(len(size)):
        text . write (str(size[i]) + "\n")
    #Write strictly 2D file
    for i in xrange(size[0]):#203
        for j in xrange(size[1]):#198
            if int(data[i][j]) == 0:
                text . write (str(-1) + "\n")
            else:
                text . write (str(int(data[i][j])) + "\n")
    text . close()

def write_sparse_array( arr, output, ndim=2 ):
    """
    Write an array to *sparse* Perseus format.
    
    arr -- numpy ndarray

    output -- name (full path) of output file

    NOTE -- 
    """
    with open( output, 'w' ) as fh:
        fh.write( str( arr.ndim )+'\n' )
        for i in xrange( arr.shape[0] ): 
            for j in xrange( arr.shape[1] ):
                if int( arr[i,j] ) != 0:
                    pos = str(i) + ' ' + str(j) + ' '
                    fh.write (pos + str(int( arr[i,j] )) + "\n")
        
    

def write_sparse_file ( fname, output=None ):
    """
    fname -- filename without extension. Extension must be .txt or .npy.

    output -- name of output file. If None, '_pers' is apended to fname, yielding <fname>_pers.txt.
    
        - Write sparse See above arguments
        """
    #dimension array
    numDim = 2
    #Open file
    if output == None:
        output = fname + '_pers'
    text = open (output + ".txt", "w")
    #Write number of dimensions
    text . write (str(numDim) + "\n")

    # now open the array file
    if fname.endswith( 'npy' ):
        data = numpy.load(fname)
    else:
        data = numpy.loadtxt(fname)
    size = data.shape
    for i in range(0,size[0]):
        for j in range(0,size[1]):
            if int(data[i][j]) != 0:
                pos = str(i) + ' ' + str(j) + ' '
                text . write (pos + str(int(data[i][j])) + "\n")
    text . close()


def write_double_file ( fname, output ):
    """
        - Write contents of single frame to output file
        e.g. new cell name on simplex:
        '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/new_110125-concatenated-ASCII_324.npy'
        e.g. old cell name on simplex:
        '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_100125/old_100125-concatenated-ASCII_1525.npy'
        Perseus input requires following format:
        #Dimension of cubical grid
        #Number of cubes in first dimension
        #Number of cubes in second dimension
        #Number of cubes in   ...     dimension
        #Birth time of cube anchored at (0,0,...,0)
        #Birth time of cube anchored at (0,0,...,1)
        #Note: Output is file name
        """
    #Load npy file
    data = numpy.load(fname)
    size = data.shape
    dim = len ( size )
    
    #Write to output
    #Open file
    text = open (output + ".txt", "w")
    #Write number of dimensions
    text . write (str(dim) + "\n")
    #Write resolutions
    for i in size:
        text . write (str(i) + "\n")
    #Write strictly 2D file
    for i in range(0,size[0]):
        for j in range(0,size[1]):
            if data[i][j] == 0:
                text . write (str(-1) + "\n")
            else:
                text . write (str(data[i][j]) + "\n")
    text . close()

def write_distance_mat( data, output=None, g=0.1, stepsize=0.2,
                        nsteps=10, dcap=None, metric='euclidean' ):
    """
    For use with

    perseus distmat <path to distance matrix file> <output string>
    """
    # load from file if data is not an array of points
    if not hasattr( data, '__index__' ):
        # .npy or .txt
        try:
            data = numpy.load( data )
        except IOError:
            data = numpy.loadtxt( data )
    dist = distance.pdist( data, metric=metric )
    if output:
        dist = distance.squareform( dist )
        space = " "
        # num points x dimension
        nx, ny = data.shape
        if not dcap:
            dcap = ny
        if not output.endswith( 'txt' ):
            output += '.txt'
        with open( output, 'w' ) as fh:
            # nx x nx distance matrix
            fh.write( str( nx )+'\n' )
            # initial threshold, step size, num steps, dimension cap==ny
            params = [ str( g ), str( stepsize ), str( nsteps ), str( dcap ) ]
            params = space.join( params )
            params += '\n'
            fh.write( params )
            for row in dist:
                r = [ str( x ) for x in row ]
                r = space.join( r )
                r += '\n'
                fh.write( r )
        return dist
    else:
        return distance.squareform( dist )

def plot_frame( frame ):
    """
        Use pylab.imshow() to plot a frame (npy array). 
        Note: The boundary should be removed, thus there will not be a halo
        """
    data = numpy.load(frame)
    fig = plt.figure()
    plt.title("RBC frame")
    plt.imshow(data)
    plt.colorbar()
    plt.show()

def plot_block( files, m, n ):
    """
        - files is directory for cell
        - m starting index, n ending index
        Use pylab.imshow() to plot a frame (npy array). 
        Note: The boundary should be removed, thus there will not be a halo
        """
    #Get all frames
    if os.path.isdir (files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'Error - input is not a directory'
        return
    #Sort frames
    frames . sort(key=natural_key)
    stack = []
    #load npy arrays
    for ind in xrange(len(frames)):
        if ind >= m and ind <= n:
            stack.append(numpy.load(files + frames[ind]))
    #Create 3d array
    d = numpy.dstack(stack)
    d = numpy.rollaxis(d,-1)
    fig = plt.figure()
    plt.title("RBC Stack")
    plt.imshow(d[1])
    plt.colorbar()
    plt.show()

def write_All_Cells ( dir, output):
    """
        Collects all cells in same directory (say '../New/frames') and writes perseus output (currently writes sparse output)
        Note: dir, output need to end with '/'
        Ex dir: '/data/rbc-Diff-Frames/New/frames/10-diff/'
        ex output: '/home/kellys/Persistence/Perseus/Formatted_Cells/frames/New/diff-10/'
    """
    cellFolders = []
    cellOutput = []
    dlist = os.listdir(dir)
    for f in dlist:
        if os.path.isdir(dir+f) and (f.startswith('new') or f.startswith('old')):
            cellFolders.append(dir+f+'/')
            cellOutput.append(output+f+'/')
    while cellFolders: #if cellF not empty
        write_Cell (cellFolders.pop(0), cellOutput.pop(0) )#call write Cell

def write_Cell ( files, output ):
    """
        - Using write_file for cell's stack of frames
        - output is location of desired output directory
        e.g. files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/'
        output = '/home/kellys/Dropbox/rbc_shared/'
        - warning: careful use of output (probably not in Dropbox)
        """
    #grab files, skip directories
    if os.path.isdir (files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'error'
    #sort frames
    frames.sort(key=natural_key)
    #for all frames call write_file
    for frame in frames:
        write_sparse_file ( fdir+frame, output+frame.rstrip('.npy'))

def write_Comp_Cell ( compArray, output ):
    """
        - WRITE COMPLIMENT CELL
        - Using write_file for cell's stack of frames
        - output is location of desired output directory
        e.g. files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/'
        output = '/home/kellys/Dropbox/rbc_shared/'
        - warning: careful use of output (probably not in Dropbox)
        """
    data = numpy.load(compArray)
    #dimension array
    numDim = 2
    #Write to output
    for ind in xrange(data.shape[0]):
        #Open file
        text = open (output +'_' + str(ind) + ".txt", "w")
        #Write number of dimensions
        text . write (str(numDim) + "\n")
        #Write strictly 2D file
        size = data[ind].shape
        for i in range(0,size[0]):
            for j in range(0,size[1]):
                if int(data[ind][i][j]) != 0:
                    pos = str(i) + ' ' + str(j) + ' '
                    text . write (pos + str(int(data[ind][i][j])) + "\n")
        text . close()
        
def write_Cell_Block ( cAbbr ):
    """
        - WRITE 3D NUMPY FILE INTO PERSEUS FORMAT
        - NOTE: cnames is a dictionary of paths to numpy files
                fnames is a dictionary of paths to output file (perseus format)
        - Using write_file for cell's stack of frames
        - output is location of desired output directory
        e.g. files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/'
        output = '/home/kellys/Dropbox/rbc_shared/'
        - warning: careful use of output (probably not in Dropbox)
        """
    cnames = pkl.load(open('moreCellNames.pkl','r'))
    fnames = pkl.load(open('moreFileNames.pkl','r'))
    data = numpy.load(cnames[cAbbr])
    output = fnames[cAbbr]
    #dimension array
    numDim = 2
    #Write to output
    for ind in xrange(data.shape[-1]):
        #Open file
        name = output.split('/')[-2] + '-concatenated-ASCII_'+str(ind)
        text = open (output + name + ".txt", "w")
        #Write number of dimensions
        text . write (str(numDim) + "\n")
        #Write strictly 2D file
        size = data[:,:,ind].shape
        for i in range(0,size[0]):
            for j in range(0,size[1]):
                if int(data[i][j][ind]) != 0:
                    pos = str(i) + ' ' + str(j) + ' '
                    text . write (pos + str(int(data[i][j][ind])) + "\n")
        text . close()
    
def write_block_dense ( files, output, m, n ):
    """
        - Write contents of 3D stack of frames to output
        required by perseus
        e.g. files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/'
        - DENSE CUBICAL TOPLEX FORMAT
        - m is index of frame to start from (index start at 0)
        - n is index of frame to end at
        """
    #grab files, skip directories
    if os.path.isdir (files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'Error - First argument not a directory!!'
        return
    #Get necessary dimensions
    numDim = 3
    dataDims = []
    dataDims . append(n-m+1)
    imgSize = numpy.load(files + frames[0]).shape
    for size in imgSize:
        dataDims . append (size)
    #Error handling for block bounds
    if n > m:
        if m < 0:
            print 'Starting index out of bounds'
            return
        if n > (len(frames)-1):
            print 'End index out of bounds'
            return
    #Sort frames
    frames.sort(key=natural_key)
    #Open file ( and write dimensions )
    text = open (output + ".txt", "w")
    text . write (str(numDim) + "\n")
    for i in xrange(len(dataDims)):
        #text . write (str(dataDims[i]) + "\n") #Iterate forward
        text . write (str(dataDims[-1-i]) + "\n")#Iterate backward
    #Write data
    for ind, frame in zip(xrange(len(frames)),frames):
        if ind >= m and ind <= n: #if within block
            data = numpy.load(files+frame)
            for i in xrange(imgSize[0]):
                for j in xrange(imgSize[1]):
                    if int(data[i][j]) == 0:
                        text . write (str(-1) + "\n")
                    else:
                        text . write (str(int(data[i][j])) + "\n")
    text . close()


def write_sparse ( files, output, m ):
    """
        - Write contents of frames to output
        required by perseus
        e.g. files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/'
        - Currently in beta
        - m is index of frame to start from (index start at 0)
        - n is index of frame to end at
        """
    #grab files, skip directories
    if os.path.isdir (files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'Error - First argument not a directory!!'
        return
    #dimension array
    numDim = 2
    nFrames = len(frames)
    dim = []
    #Error handling for bounds
    if m < 0 or m > (nFrames-1):
        print 'index out of bounds'
        return
    #append number of frames, Perseus doesn't use this
    #dim.append( nFrames )
    size = numpy.load(files + frames[0]).shape
    #append sizes of frame, Perseus doesn't use this
    #dim.append(size[0])
    #dim.append(size[1])
    #Sort frames
    frames.sort(key=natural_key)
    #Write to output
    #Open file
    text = open (output + ".txt", "w")
    #Write number of dimensions
    text . write (str(numDim) + "\n")
    #Write strictly 2D file
    data = numpy.load(files + frames[m])
    for i in xrange(size[0]):
        for j in xrange(size[1]):
            if int(data[i][j]) != 0:
                pos = str(i) + ' ' + str(j) + ' '
                text . write (pos + str(int(data[i][j])) + "\n")
    text . close()

def write_block_sparse ( files, output, m, n ):
    """
        - Write contents of 3D stack of frames to output
        required by perseus
        e.g. files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/'
        - Currently in beta
        - m is index of frame to start from (index start at 0)
        - n is index of frame to end at
        """
    #grab files, skip directories
    if os.path.isdir (files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'Error - First argument not a directory!!'
        return
    #dimension array
    numDim = 3
    nFrames = len(frames)
    dim = []
    #Error handling for block bounds
    if n > m:
        if m < 0:
            print 'Starting index out of bounds'
            return
        if n > (nFrames-1):
            print 'End index out of bounds'
            return
    #dim.append( nFrames )
    size = numpy.load(files + frames[0]).shape
    #Sort frames
    frames.sort(key=natural_key)
    #Write to output
    #Open file
    text = open (output + ".txt", "w")
    #Write number of dimensions
    text . write (str(numDim) + "\n")
    #Write strictly 2D file
    for ind, frame in zip(xrange(len(frames)),frames):
        if ind >= m and ind <= n: #if within block
            data = numpy.load(files + frame)
            for i in xrange(size[0]):
                for j in xrange(size[1]):
                    if int(data[i][j]) != 0:
                        pos =  str(ind) + ' ' + str(i) + ' ' + str(j) + ' ' 
                        text . write ( pos + str(int(data[i][j])) + "\n")
    text . close()

def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

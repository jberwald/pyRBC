import numpy
import time
import csv
import os


def cell2npy( fname, skip=0 ):
    """
    Sometimes a cell image is stored in a DOS file format, with
    approximately 200 lines. This reads the file and converts the
    values to ints before stacking them in a numpy array. The first
    two lines are the dimension of the cell matrix. The third line is
    a mysterious 'I' character. The fourth entry is a blank line. We
    skip these and fill the numpy arrays from the fifth line (i=4 in
    zero-based indexing).

    fname -- full path to the file
    """
    # this convert from DOS newline format ('U'==universal newline)
    reader = csv.reader( open( fname, 'rU' ), delimiter='\t' )
    lines = list( reader )
    # convert strings to ints
    rows = [ [ int( val ) for val in row[:-1] ] for row in lines[skip:] ]
    try:
        return numpy.vstack( rows[skip:] )
    except ValueError:
        print "problem with rows, here are the first 10:", rows[:10]
                       
    
def stack_cell( cell, fdir ):
    """
    fdir -- directory containing the Complement cell image frames.

    Return a 3D numpy array, with the third axis corrseponding to
    time.
    """
    if not fdir.endswith( '/' ): fdir+='/'
    dlist = os.listdir( fdir + cell + '/' )
    dlist.sort()
    all_frames = []
    print "converting frames..."
    tstart = time.time()
    for frame in dlist:
        # this skips any hidden files... a little error handling at
        # least :}
        if frame.startswith('.'): continue
        all_frames.append( cell2npy( fdir + cell + '/' + frame ) )
    print "It took ", time.time() - tstart, "seconds for converting", \
        len( dlist ), "frames."
    numpy.save( fdir + cell, all_frames )
    return numpy.array( all_frames, dtype=numpy.int )

if __name__ == "__main__":

    celldir = '/data/jberwald/wyss/data/rbc/cells/'
    cells = [ 'new_3', 'new_4', 'new_5', 'new_6', 'new_7', 'new_8', 'new_9',
              'old_10', 'old_2', 'old_5', 'old_6', 'old_7', 'old_8', 'old_9' ]

    for cell in cells:
        print "Lumping frames for cell", cell, "into a 3D array..."
        stack_cell( cell, celldir )

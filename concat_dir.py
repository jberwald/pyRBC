import numpy as np
import os, time, shutil
import cPickle as pkl
import fft_image as F


def load_pkl( frame ):
    with open( frame ) as fh:
        f = pkl.load( fh )
    return f

def concat_files( fdir ):
    """
    This assume we are on a Linux box (i.e. Sciclone). Concatenate all
    files in a threshold directory for easier analysis.
    """
    tstart = time.time()
    if not fdir.endswith( '/' ): fdir += '/'
    frames = os.listdir( fdir )
    # if os.uname()[0] == 'Linux':
    #     frames = []
    #     for f in dlist:
    #         if f.endswith('npy') and not os.path.isdir( fdir+f ):
    #             frames.append( f )
    #         else:
    #             frames = [ f for f in dlist if not os.path.isdir( fdir+f ) ]

    if frames[0].endswith( 'pkl' ):
        load_frame = load_pkl
    else:
        load_frame = np.loadtxt
    # dict for ordered hash map
    fd = {}
    # order frames numerically
    for f in frames:
        sf = f.split( '_' )
        try:
            fd[ int( sf[-2] ) ] = load_frame( fdir + f )
        # some cubical files are corrupted, probably due to thresholding issue. Must look into this. 
        except ValueError:
            print f
            continue
    savename = fdir + 'all_frames.pkl'
    with open( savename, 'w' ) as fh:
        pkl.dump( fd, fh )
    total_time = time.time() - tstart
    print "Saved files to", savename, " in ", total_time, " seconds."


def copy_concat_file( cell, mode,
                       new_home='/sciclone/data10/jberwald/wyss/concat_files/',
                       old_home='/sciclone/data10/jberwald/wyss/data/Cells_Jesse/'):
    """
    Copy the concatenated fft files (at a given mode) to a single
    directory. Rename with cell and mode appended in the process.

    Files copied to <new_home>/cell/r<mode>_betti/all_pickle.pkl
    """
    if not cell.endswith( '/' ): cell += '/'
    # make a new directory to copy file to 
    F.make_dir( new_home + cell )
    mode_str = F.mode2str( mode )
    r_mode = 'r' + mode_str + '_betti/'
    new_fpath = new_home + cell + r_mode 
    F.make_dir( new_fpath )

    # now copy from old_home -> new_home
    if 'old' in cell:
        cell_type = 'Old/'
    else:
        cell_type = 'New/'
    fft_path = 'fft_frames/normed_frames/'
    fpath = old_home + cell_type + 'frames/' + cell + fft_path +\
        r_mode + 'all_frames.pkl'
    shutil.copy( fpath, new_fpath ) 

if __name__ == "__main__":

    modes = np.linspace( 0, 1, 21 )

    new_fdir = '/sciclone/home04/jberwald/data10/wyss/data/Cells_Jesse/New/frames/'
    old_fdir = '/sciclone/home04/jberwald/data10/wyss/data/Cells_Jesse/Old/frames/'
    new_cells =[  'new_110125/' ]#, 'new_130125/', 'new_140125/', 'new_40125/', 'new_50125/' ] #   'new_110125/',
    old_cells = [  'old_100125/',  'old_120125/', 'old_50125/', 'old_90125/' ] #,

    subdir = 'fft_frames/normed_frames/'

    print "concatenating cell data..."
    for cell in new_cells:
        print "cells", cell
        #concat_files( new_fdir _ cell )
        for mode in modes:
            str_mode = F.mode2str( mode )
            #thresh_dir = 'r' + str_mode + '/'
            betti_dir = 'r' + str_mode + '_betti/'
            # first the threshold files
            # cell_dir = new_fdir + cell + subdir + thresh_dir
            # concat_files( cell_dir )
            # now the betti files
            cell_dir = new_fdir + cell + subdir + betti_dir
            concat_files( cell_dir )
        
    
    # for cell in old_cells:
    #     print "cells", cell
    #     #concat_files( new_fdir _ cell )
    #     for mode in modes:
    #         #thresh_dir = 'r' + str_mode + '/'
    #         str_mode = F.mode2str( mode )
    #         betti_dir = 'r' + str_mode + '_betti/'
    #         # first the threshold files
    #         # cell_dir = old_fdir + cell + subdir + thresh_dir
    #         # concat_files( cell_dir )
    #         # now the betti files
    #         cell_dir = old_fdir + cell + subdir + betti_dir
    #         concat_files( cell_dir )
        

    # for cell in old_cells:
    #     print "cells", cell
    #     #concat_files( new_fdir _ cell )
    #     for mode in modes:
    #         copy_concat_file( cell, mode )

    print "copying concatentated files..."
    for cell in new_cells:
        print "cells", cell
        #concat_files( new_fdir _ cell )
        for mode in modes:
            copy_concat_file( cell, mode )

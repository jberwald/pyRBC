"""
    
Module to run perseus
created by kel 5/23/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""

import numpy
import matplotlib.pylab as plt
import re
import os
import rbc_npy2Perseus as R
from multiprocessing import Pool
from subprocess import call
import pickle as pkl

def perseus ( fname, output, type='scubtop' ):
    """
        TYPE is input to perseus ex, cubtop, scubtop, etc
        run perseus
    """
    call(['/usr/bin/perseus', type, fname, output])

def run_perseus_All ( type, dir, output):
    """
        Collects all cells in same directory (say '../New') and runs perseus
        Note: dir, output need to end with '/'
    """
    cellFolders = []
    cellOutput = []
    dlist = os.listdir(dir)
    for f in dlist:
        if os.path.isdir(dir+f):
            cellFolders.append(dir+f+'/')
            cellOutput.append(output+f+'/')#no +f if want in ..-all directory
    while cellFolders:
        run_perseus (type, cellFolders.pop(0), cellOutput.pop(0) )#perseus for frames of cell

def run_perseus_More ( type='scubtop' ):
    """
        Collects all cells in same directory (say '../New') and runs perseus
        Note: dir, output need to end with '/'
    """
    fnames = pkl.load(open('moreFileNames.pkl','r'))
    outnames = pkl.load(open('morePersFiles.pkl','r'))
    
    for cAbbr in fnames:
        run_perseus (type, fnames[cAbbr], outnames[cAbbr] )#perseus for frames of cell

def run_perseus_All_parallel ( dir, output):
    """
        Parallel implementation
        """
    pool = Pool(processes=7)
    cellFolders = []
    cellOutput = []
    dlist = os.listdir(dir)
    for f in dlist:
        if os.path.isdir(dir+f):
            cellFolders.append(dir+f+'/')
            cellOutput.append(output+'/')#no +f if want in ..-all directory
    pool.map(run_perseus,zip(cellFolders,cellOutput))#perseus for frames of cell

def run_perseus_list ( type, folder, file, list, output):
    """
        folder - directory to cells, end in '/'
        file - category of cell ex: new_110125
        list of frames to calculate
        output directory, end in '/'
    """
    interm = '-concatenated-ASCII_'
    for l in list:
        dir = folder + file + '/'
        frame = dir + file + interm + str(l) + '.txt'
        perseus ( type, frame, output + file + '_' + str(l))

def run_perseus ( type, files, output ):
    """
        - Using perseus to run on single frame
        - output is location of desired output directory
        e.g. files = '/home/kellys/Persistence/Perseus/Output/frames/Old/diff-10/'
        output = '/home/kellys/Persistence/Perseus/Output'
        - warning: careful use of output (probably not in Dropbox)
        """
    #grab files, skip directories
    if os.path.isdir (files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('.txt') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'error'
    #sort frames
    frames.sort(key=R.natural_key)
    #for all frames call write_file
    for frame in frames:
        print fdir+frame
        perseus ( type, fdir+frame, output+frame.rstrip('.txt'))
        

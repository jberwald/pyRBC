"""
    
Module to computing clusters of persistence diagrams, also various functions for related tasks

    e.g. files on simplex:
        '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/'    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""

import hcluster
import numpy
import scipy.cluster.hierarchy as hier
import matplotlib.pylab as plt
import re
import os
from subprocess import call
import scipy.io
import itertools


def dmat_parallel (cells, files, frameList, output, type, b_num):
    """
        ---UNDER CONSTRUCTION---
        Compute distance matrix in parallel
        arguments described in function below
        """
    nFrames = len(cells)*len(frameList)
    data = numpy.zeros((nFrames,nFrames), dtype=numpy.uint16)
    cellFrames = []
    for i in xrange(len(cells)):
        dlist = os.listdir(cells[i])
        #for each frame
        for f in dlist:
            #if ends with correct betti number
            if f.endswith(str(b_num)+'.txt'):
                #if the index of array (ex .._99_..) is in frameList
                if numpy.where(frameList == int(f.split('_')[2]))[0].size > 0:
                    trueVal = 0
                    # check whether desired category (ex new_110125)
                    for file in files:
                        if f.split('-')[0] == file.split('/')[-2]:
                            trueVal = 1
                    if trueVal == 1: #if true append
                        cellFrames.append(f)
    #sort the frames
    cellFrames.sort(key=natural_key)
    cellList = [cellFrames, cellFrames]
    cartProd = list(itertools.product(*cellList))
    #multiprocessing here
    dist = persistence_distance(type,persOut,persOut2,fname,fname2)
    data[i][j] = dist
    data[j][i] = dist
    saveMatrix(data, output)

def dmat (cells, files, frameList, type, b_num, save='NO', output='NIL'):
    """
        Simplified distance matrix creation
        Argument description:
        cells -  array of paths to cells for each cell, ex:
        '/data/PerseusOutput/..'
        files - array of paths to corresponding data files (needed for getting max height values)
        '/data/jberwald/wyss/...'
        frameList - numpy array of desired frame indices (start at 0)
        """
    nFrames = len(files)*len(frameList)
    data = numpy.zeros((nFrames,nFrames), dtype=numpy.uint16)
    #list for each frame's perseus output
    cellFrames = []
    #dictionary for frame locations (to get max height value)
    cellDict = {}
    #for each cell or 'category' specified
    for i in xrange(len(cells)):
        dlist = os.listdir(cells[i]) #get all cells in directory
        #for each frame
        for f in dlist:
            #if ends with correct betti number
            if f.endswith(str(b_num)+'.txt'):
                #if the index of array (.._99_..) is in frameList
                if numpy.where(frameList == int(f.split('_')[2]))[0].size > 0:
                    # check whether desired category (ex new_110125), if so get the index
                    frameNames = [files[x].split('/')[-2] for x in xrange(len(files))]
                    index = [j for j, x in enumerate(frameNames) if x == f.split('-')[0]]
                    if index:
                        #put onto cellFrames, create mapping from cellFrame to its associated location
                        cellPath = cells[i] + f
                        fname = files[index[0]] + f[:-6] + '.npy'
                        cellFrames.append(cellPath)
                        cellDict[cellPath] = fname
    cellFrames.sort(key=natural_key)
    for f, g in itertools.product(enumerate(cellFrames),enumerate(cellFrames)):
        (f_ind, f) = f
        (g_ind, g) = g
        if g_ind >= f_ind:  #exploiting symmetry
            dist = persistence_distance(type,f,g,cellDict[f],cellDict[g])
            data[f_ind][g_ind] = dist
            data[g_ind][f_ind] = dist
    if save == 'YES':
        saveMatrix(data,output)
    return data

def dvec ( cell, file, frameList, type, b_num, save='NO', output='NIL'):
    """
        For fixed frame, calculate distance vector between other frames (known as v_k)
        Arguments:
        cell is directory path to cell location
        file is location of specific specific frame for v_k
        ex: '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/new_110125-concatenated-ASCII_99.npy'
        frameList is list of indices to calculate against
        output is location (name inclusive) for saving
        type is metric 'bd' or 'wass', b_num is betti number
    """
    #Get length of vector
    #check if index '..._index.npy' is in frameList
    index = file.rstrip('.npy').split('_')[-1]
    if index not in frameList:  #add fixed index if not in list of frames
        frameList = numpy.append(index)
    nFrames = len(frameList)
    #sort indices
    frameList.sort()
    #create data vector
    data = numpy.zeros(nFrames)
    cellFrames = [] #list of paths to the frames in frameList
    cellDict = {}   #mapping to npy arrays for each frame
    indexCellPath = ''  #cell path for specified fixed file
    #strips 'new_11....npy' from file to get general path
    filePath = file.rstrip(file.split('/')[-1])
    dlist = os.listdir(cells)
    for f in dlist:
        if f.endswith(str(b_num)+'.txt'):#correct betti number
            #get index of frame, i.e. '..._99.npy'
            frameNum = int(f.split('_')[2])
            if frameNum in frameList:   #correct frame index
                #correct type, ex. new_110125
                if f.split('-')[0] == file.split('/')[-2]:
                    cellPath = cell + f
                    fname = filePath + f[:-6] + '.npy'
                    cellFrames.append(cellPath)
                    cellDict[cellPath] = fname
                    #if specified fixed file
                    if frameNum == index:
                        indexCellPath = cellPath
    cellFrames.sort(key=natural_key)
    for g_ind, g in enumerate(cellFrames):
        # Compute distance
        data[g_ind] = persistence_distance(type,indexCellPath, g, cellDict[indexCellPath], cellDict[g])
    if save=='YES':
        saveMatrix(data,output)
    return data

def dlag_vec (cell, files, lag, type, b_num, tempOut='/data/tOut.txt', save='NO', output='NIL'):
    """
        Compute lag vector
        Arguments are similar as above
        ex cell = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/new_110125/'
        BUT hanldes the case '/data/PerseusOutput/original/2d_sparse/New/' (i.e. All in one file)
        files = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/'
        type = 'bd' or 'wass'
        Except for lag, which is the difference variable for computing v_lag ex. lag = 25
    """
    cellFrames = []
    cellDict = {}
    dlist = os.listdir(cell)
    for f in dlist:
        if f.endswith(str(b_num)+'.txt'):    #correct betti num
            if f.split('-')[0] == files.split('/')[-2]:   #correct file type
                cellPath = cell + f
                fname = files + f[:-6] + '.npy'
                cellFrames.append(cellPath)
                cellDict[cellPath] = fname
    cellFrames.sort(key=natural_key)
    data = numpy.zeros(len(cellFrames) - (lag-1))#allocate length based on lag
    cellStack = []
    for ind, g in enumerate(cellFrames):
        if ind % lag == 0: #TESTING PURPOSES --- GETTING MUCH SMALLER VECTOR
            if ind >= lag:
                f = cellStack.pop(0)
                data[ind-lag] = persistence_distance(type,f,g,cellDict[f],cellDict[g],tempOut)
            cellStack.append(g)
    if save=='YES':
        saveMatrix(data,output)
    return data
                               
def persistence_distance (type, persFile1, persFile2, fname1, fname2,tempOut, p=2):
    """
        Compute persistence distance ('bd' or 'wass')
        Call modified metric code in kel_metric
        Preemptive get max heights from .npy for infinite persistence
        p for pth Wasserstein distance (default of 2)
    """
    maxl = get_Max(fname1)
    maxl2 = get_Max(fname2)
    print '\n...computing ' + type + ' distance..'
    print persFile1
    print persFile2
    if type.startswith('b'):
        call(['/home/kellys/Dropbox/WM/KellySpendlove/metric/kel_metric/bottleneck/bd', persFile1, persFile2, str(maxl), str(maxl2), tempOut])
    else:
        call(['/home/kellys/Dropbox/WM/KellySpendlove/metric/kel_metric/wasserstein/src/wass', persFile1, persFile2, str(p), str(maxl), str(maxl2), tempOut])
    with open(tempOut, 'r') as f:
        d = float(f.read())
    f.closed
    os.remove(tempOut)
    return d

def get_Norm ( data, cast ):
    """
        Compute norm, cast is norm desired
        Ex. cast = numpy.inf (max norm), cast = 1 (1-norm), cast = 2 (2-norm)
        cast = 'fro' (frobenius)
    """
    return numpy.linalg.norm(data, cast)

def get_Max ( fname, add=0 ):
    """
        Send in add=k to return max height value + k
        For example, setting add=1 works when H_1 
    """
    return numpy.load(fname).max()+add

def saveMatrix ( data, output ):
    scipy.io.savemat(output + '.mat', mdict={'arr':data})

def normalize ( frame ):
    """ normalize cell """
    max = frame.max()
    d = frame.astype(float) / max
    return d

def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

#Deprecated Functions 6/12/2012

def create_dist_matrix (cells, files, frameList, output, type, b_num):
    """
        DEPRECATED (USE DMAT).
        POTENTIAL BUG.
        cells -  array of paths to cells for each cell, ex:
        '/data/PerseusOutput/..'
        files - array of paths to corresponding data files
        '/data/jberwald/wyss/...'
        frameList - list of desired frame indices (start at 0)
        --------DEPRECATED-----
        """
    nFrames = len(cells)*len(frameList)
    data = numpy.zeros((nFrames,nFrames), dtype=numpy.uint16)
    for i in xrange(len(cells)):
        dlist = os.listdir(cells[i])
        cellFrames = []
        #for each frame
        for f in dlist:
            #if ends with correct betti number
            if f.endswith(str(b_num)+'.txt'):
                #if the index of array (.._99_..) is in frameList
                if numpy.where(frameList == int(f.split('_')[2]))[0].size > 0:
                    trueVal = 0
                    # check whether desired category (ex new_110125)
                    for file in files:
                        if f.split('-')[0] == file.split('/')[-2]:
                            trueVal = 1
                    if trueVal == 1: #if true append
                        cellFrames.append(f)
        cellFrames.sort(key=natural_key)
        for j in xrange(len(cells)):
            if j >= i:
                dlist2 = os.listdir(cells[j])
                cellFrames2 = []
                for f in dlist:
                    if f.endswith(str(b_num)+'.txt'):
                        if numpy.where(frameList == int(f.split('_')[2]))[0].size > 0:
                            trueVal = 0
                            for file in files:
                                if f.split('-')[0] == file.split('/')[-2]:
                                    trueVal = 1
                            if trueVal == 1:
                                cellFrames2.append(f)
                cellFrames2.sort(key=natural_key)
                #get index and frame
                for ind, frame in zip(xrange(len(cellFrames)),cellFrames):
                    persOut = cells[i] + frame
                    fname = files[i] + frame[:-6] + '.npy'
                    #get index and frame of associated cells
                    for ind2, frame2 in zip(xrange(len(cellFrames2)),cellFrames2):
                        persOut2 = cells[j] + frame2
                        fname2 = files[j] + frame2[:-6] + '.npy'
                        dist = persistence_distance(type,persOut,persOut2,fname,fname2)
                        data[ind*i][ind2*j] = dist
                        data[ind2*j][ind*i] = dist
    saveMatrix(data, output)


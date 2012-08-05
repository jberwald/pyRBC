"""
    
Module to convert RBC frames (npy files) to input required by CHoMP
Module for running CHomP
    - 
    created by kel 6/23/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
    
    
"""

import numpy
import matplotlib.pylab as plt
import re
import subprocess, os
import pickle as pkl

def npy3d2chomp ( arr, output, threshold ):
    """
        Convert numpy array to CHomP tuple format
    """
    aMean = arr[arr>0].mean()
    tHeight = threshold*aMean
    tupleList = []
    for i in arr.shape[0]:
        for j in arr.shape[1]:
            for k in arr.shape[2]:
                if arr[i][j][k] <= tHeight:
                    tupleList.append((i,j,k))
    tupleArr = numpy.array(tupleList)
    array2chomp(tupleArr,output)

def npy2chomp ( arr, k, output, threshold ):
    """
        Convert numpy array to CHomP tuple format
        """
    aMean = arr[arr>0].mean()
    tHeight = threshold*aMean
    tupleList = []
    for i in xrange(arr.shape[0]):
        for j in xrange(arr.shape[1]):
            if arr[i][j] <= tHeight and not int(arr[i][j])==0:
                tupleList.append((k,i,j))
    tupleArr = numpy.array(tupleList)
    array2chomp(tupleArr,output)

def array2chomp( arr, savename ):
    """
        Convert an array to chomp format, ( , , ). Write the resulting
        column of numbers to disk.
        """
    rows = map( lambda x: str(x)+'\n', map( tuple, iter( arr ) ) ) 
    with open( savename, 'a' ) as fh:
        fh.writelines( rows )
    fh.close()

def run_chomp( fname, savename ):
    """
        Call chomp to compute the betti numbers of the image in file fname.
        
        See http://chomp.rutgers.edu
        """
    cmd = ['chomp', fname]
    try:
        with open( savename, 'w' ) as fh:
            p = subprocess.call( cmd, stdout=fh )
    except:
        print "subprocess returned with code", p

def chomp_block (folder,output,m,n, threshold):
    #grab files, skip directories
    if not folder.endswith('/'):
        folder+='/'
    if os.path.isdir (folder):
        dlist = os.listdir(folder)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(folder+f):
                frames.append(f)
    else:
        print 'Error - First argument not a directory!!'
        return
    frames.sort(key=natural_key)
    for ind, frame in zip(xrange(len(frames)),frames):
        if ind >= m and ind <= n:
            data = numpy.load(folder+frame)
            npy2chomp(data, ind, output, threshold)

def run_block(cell,m,n,thresholds):
    """
        Run chomp on 3d block
        cell is of form 'n11', 'o9', etc...
        m, n starting, ending indices
        thresholds list of thresholds
    """
        
    fnames = pkl.load(open('fileNames.pkl','r'))
    folder = fnames[cell]
    path = folder.split('Cells_Jesse')[-1]
    title = path.split('/')[-2] #ex 'new_110125/'
    subtitle = '_'+str(m)+'_'+str(n)+'_'
    cubStub = '/data/ChompData/Blocks'+path+title+subtitle
    outStub = '/data/ChompData/ChompOutput'+path+title+subtitle
    
    for t in thresholds:
        str_t = ''.join(str(t).split('.'))
        cubFile = cubStub+str_t+'.cub'
        chomp_block(fnames[cell],cubFile,m,n,t)
        savename = outStub+str_t+'.txt'
        run_chomp(cubFile,savename)

def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

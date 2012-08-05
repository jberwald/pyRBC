"""
    
    Extract betti number time series for each frame in cell (for specified threshold)
    
    Uses the output of Perseus in order to avoid having to recompute 
    the homology of various thresholds
   
created by kel 6/22/2012
    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
    
    
"""

import numpy
import matplotlib.pyplot as plt
import re
import os
import rbc_npy2Perseus as n2p

def extract_betti (persFiles, files, b_num, threshold):
    """
        Extract betti number time series from cell at for sublevel sets below threshold
        persFiles is directory of perseus output for cell
        files is directory for associated .npys
        b_num is betti number
        threshold is factor of mean desired (i.e. .5-1.5)
    """
    if not persFiles.endswith('/'):
        persFiles+='/'
    if not files.endswith('/'):
        files+='/'
    if os.path.isdir (persFiles):
        dlist = os.listdir(persFiles)
        frames = []
        for frame in dlist:
            if frame.endswith('betti.txt') and not os.path.isdir(persFiles+frame):
                frames.append(frame)
    else:
        print 'Error: Not a Directory'
        return
    frames.sort(key=n2p.natural_key)
    betti_ts = []
    betti_heights = []
    betti_list = []
    betti_list.append(betti_heights)
    betti_list.append(betti_ts)
    for ind, frame in zip(xrange(len(frames)),frames):
        with open(persFiles+frame, 'r') as f:
            s = f.read()
        f.close()
        data = numpy.load(files+frame.rstrip('_betti.txt')+'.npy')
        #dM = data[data>0].mean() #USE MEAN
        dM = numpy.median(data[data>0])
        tHeight = threshold*dM
        dstr = s.split('\n')
        dstr.pop(0)#pop '' - could use .remove('')
        for j, subs in zip(xrange(len(dstr)),dstr):
            curHeight = int(re.findall('\[[^\]]*\]',subs)[0][1:-1])
            if tHeight < curHeight:
                betti_str = dstr[j-1].split(':')[-1][1:-1]
                bHeight = int(re.findall('\[[^\]]*\]',dstr[j-1])[0][1:-1])
                betti_ts.append (betti_str.split(' ')[b_num])
                betti_heights.append ( bHeight )
                break
    if len(betti_ts) != len(frames):
        print 'ERROR IN BETTI TIME SERIES SIZE - TRY DIFFERENT THRESHOLD VALUE'
    return betti_list
                         
def plot_betti (persCells, cells, b_num, threshold):
    """
        Plot betti number time series
    """
    if len(persCells) != len(cells):
        print 'INCORRECT INPUT'
        return
    ts = []
    tHeights = []
    colorArr = ['b', 'r', 'g', 'y', 'c', 'm']
    for i in xrange(len(persCells)):
        blist = extract_betti(persCells[i],cells[i], b_num,threshold)
        tHeights.append(blist[0])
        ts.append( blist[1] )
    fig = plt.figure()
    for ind in xrange(len(ts)):
        plt.plot(ts[ind], colorArr[ind])
    plt.show()
    return tHeights, ts
        
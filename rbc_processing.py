"""
    
Module for assorted processing of RBCs and diff images
created by kel 5/23/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os

def save_npy_as_png( data, output ):
    plt.ioff()
    """
        save as npy as png
    """
    fig = plt.figure()
    plt.title('RBC ' + output)
    plt.imshow(data)
    plt.colorbar()
    plt.savefig(output)

def npy2png( folder, file, list, output ):
    plt.ioff()
    """
        folder - directory to cells
        file - category of cell ex: new_110125
        list of frames
        output directory, need end in '/'
        """
    interm = '-concatenated-ASCII_'
    for l in list:
        dir = folder + file + '/'
        frame = dir + file + interm + str(l) + '.npy'
        data = numpy.load(frame)
        fig = plt.figure()
        plt.title("RBC frame " + folder+file + ' ' + str(l))
        plt.imshow(data)
        plt.colorbar()
        plt.savefig(output + file + '_' + str(l))

def diff_image (files, savedir):
    if os.path.isdir(files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'error'
    k = 0
    n = len(frames)
    derivArr = []
    parArr = []
    frames.sort(key=natural_key)
    for frame in frames:
        k = k+1
        newFrame = numpy.load(fdir+frame)
        savename = frame.rstrip('.npy')
        if k!= 1 and k < n:
            diffFrame = numpy.abs(newFrame-prevFrame)
            numpy.save(savedir+savename + '.npy',diffFrame)
        prevFrame = newFrame

def k_diff_image (files, savedir, k):
    """
        savedir needs to be directory
        """
    if os.path.isdir(files):
        fdir = files + '/'
        dlist = os.listdir(fdir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(fdir+f):
                frames.append(f)
    else:
        print 'error'
    i = 0
    n = len(frames) - (k-1)
    frameStack = []
    frames.sort(key=natural_key)
    for frame in frames:
        i = i+1
        newFrame = numpy.load(fdir+frame)
        savename = frame.rstrip('.npy')
        if i > k:
            #diffFrame = numpy.abs(newFrame-frameStack.pop(0))#abs value
            #make some masks
            newfMask = numpy.ma.masked_less_equal(newFrame,0)
            fsMask = numpy.ma.masked_less_equal(frameStack.pop(0),0)
            diffFrame = newfMask - fsMask
            diffFrame = diffFrame + abs(diffFrame.min())#add the minimum value
            numpy.save(savedir+savename + '.npy',diffFrame.data)
        frameStack.append(newFrame)

def run_diff_All ( dir, output, k):
    """
        Collects all cells in same directory (say '../New') and runs perseus
        Note: dir, output need to end with '/'
        """
    cellFolders = []
    cellOutput = []
    dlist = os.listdir(dir)
    for f in dlist:
        if os.path.isdir(dir+f) and (f.startswith('new') or f.startswith('old')):
            cellFolders.append(dir+f+'/')
            cellOutput.append(output+f+'/')
    while cellFolders:
        k_diff_image(cellFolders.pop(0), cellOutput.pop(0), k )#perseus for frames of cell

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

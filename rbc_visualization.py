"""
    
Module to convert RBC frames (npy files) to input required by Perseus software (for computing topological persistence)
created by kel 6/5/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools as it
import re
import scipy.io
import os
import cPickle as pkl
import rbc_processing as rp

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

def plot_sublevels (file, persFile, output, interval, mask='n'):
    """
        plots sublevel sets for specified numpy file
        from the height values in associated persFile
        output is directory for output files (indv. names come from file + height)
        interval provides info on what frames to save
            (i.e. interval==2 means every other frame )
        if mask=='y' all sublevel values are set to 1
    """
    if not output.endswith('/'):
        output+='/'
    with open(persFile, 'r') as f:
        s = f.read()
    f.close()
    data = numpy.load(file)    
    heightList = [x[1:-1] for x in re.findall('\[[^\]]*\]',s)]    
    #d = numpy.where(data>int(heightList[100]))
    #data[d] = 0
    #rp . save_npy_as_png(data,output+file.split('/')[-1].rstrip('.npy'))
    for ind, h in enumerate(heightList):
        if ind % interval == 0:
            temp = data.copy()
            temp[numpy.where(temp>int(h))] = 0
            if mask == 'y':
                temp[numpy.where(temp!=0)] = 1
            outName = file.split('/')[-1].rstrip('.npy') + '_' + h
            rp . save_npy_as_png ( temp, output + outName)

def plot_sublevels_comp (file, slice, output, interval, persFile='NULL', mask='NO'):
    """
        plots sublevel sets for specified numpy file
        from the height values in associated persFile
        output is directory for output files (indv. names come from file + height)
        interval provides info on what frames to save
        (i.e. interval==2 means every other frame )
        if mask=='y' all sublevel values are set to 1
        """
    if not output.endswith('/'):
        output+='/'
    data = numpy.load(file)[slice]
    if persFile != 'NULL':
        with open(persFile, 'r') as f:
            s = f.read()
        f.close()
        heightList = [x[1:-1] for x in re.findall('\[[^\]]*\]',s)]    
    else:
        tempMin = data[numpy.where(data>data.min())].min()
        avg = (data.max() - tempMin) / 50
        heightList = numpy.arange(tempMin,data.max(), avg)
    #d = numpy.where(data>int(heightList[100]))
    #data[d] = 0
    #rp . save_npy_as_png(data,output+file.split('/')[-1].rstrip('.npy'))
    for ind, h in enumerate(heightList):
        if ind % interval == 0:
            temp = data.copy()
            temp[numpy.where(temp>int(h))] = 0
            if mask == 'YES':
                temp[numpy.where(temp!=0)] = 1
            outName = file.split('/')[-1].rstrip('.npy') + '_' + str(slice) +'_' + str(h)
            rp . save_npy_as_png ( temp, output + outName +'.png')

def plot_interlevels (file, persFile, output, interval, r, mask='NO'):
    """
        plots interlevels sets for specified numpy file
        from the height values in associated persFile (USE .BETTI FILE)
        
        output is directory for output files (indv. names come from file + height)
        interval provides info on what frames to save
        (i.e. interval==2 means every other frame )
        
        r perburation to height function (use half of desired window)
        Use half interval to capture every height value (if all height values present)
        
        if mask=='YES' all sublevel values are set to 1
        """
    if not output.endswith('/'):
        output+='/'
    with open(persFile, 'r') as f:
        s = f.read()
    f.close()
    data = numpy.load(file)
    heightList = [x[1:-1] for x in re.findall('\[[^\]]*\]',s)]    
    #d = numpy.where(data>int(heightList[100]))
    #data[d] = 0
    #rp . save_npy_as_png(data,output+file.split('/')[-1].rstrip('.npy'))
    for ind, h in enumerate(heightList):
        if ind % interval == 0:
            temp = data.copy()
            lower = int(h) - r
            upper = int(h) + r
            temp[(temp>upper)|(temp<lower)] = 0
            if mask == 'YES':
                temp[numpy.where(temp!=0)] = 1
            outName = file.split('/')[-1].rstrip('.npy') + '_'+h+'-'+str(lower)+'-'+str(upper)
            rp . save_npy_as_png ( temp, output + outName + '.png')

def create_3d_sublevels( folder, file, frameList, heightList, output, ext='npy' ):
    """
        This function creates .mats to view with Matlab or .npy
        - folder is directory to cells, ends in '/'
        - file is category of cell, ex. new_110125
        - frameList is numbers of frames to calculate on
        - heightList is list of heights for sublevel sets
        - output is output file for .mat files
        - ext is extension '.npy' or 'm'
        """
    #correcting data
    if folder.endswith('/') == False:
        folder = folder + '/'
    #processing variables
    interm = '-concatenated-ASCII_'
    # for each height
    for h in heightList: 
        #stack frames
        stack = []
        for f in frameList:
            #frame processing
            dir = folder + file + '/'
            frame = dir + file + interm + str(f) + '.npy'
            data = numpy.load(frame)
            #find superlevel set and set to 0
            d = numpy.where(data>h)
            data[d] = 0
            #place on stack
            stack . append(data)
        sublevel = numpy.dstack(stack)
        #setting dimensions
        sublevel = numpy.rollaxis(sublevel,-1)
        #save for visualization in sage or matlab
        if ext == 'npy':
            numpy.save(output + '_height_' + str(h), sublevel)
        elif ext == 'pkl':
            with open( output + '_height_' + str(h), 'w') as fh:
                pkl.dump( sublevel, fh )
        else:
            scipy.io.savemat(output + '.mat', mdict={'arr':sublevel})

def plot_frame_3d ( folder, file, list, output):
    """
        Use pylab.imshow() to plot a wire frame of npy array,
        where height > 0 in 3 dimensions.
        """
    interm = '-concatenated-ASCII_'
    if not file.endswith('/'):
        file+=('/')
    for l in list:
        fig = plt.figure()
        data = numpy.load(folder+file+interm+str(l)+'.npy')
        d = numpy.where(data>0)
        xlist = d[0]
        ylist = d[1]
        z = data[d]
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.title('RBC Frame 3d ' + frame)
        ax.plot(xlist,ylist,z)
        plt.savefig(output + file + '_' + str(l))
    
def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

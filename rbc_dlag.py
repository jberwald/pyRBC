"""
    
Module to computing lag vectors
    e.g. files on simplex:
        '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/'    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""
import numpy
import scipy.cluster.hierarchy as hier
import re
import os
from subprocess import call
import scipy.io
import itertools
import rbc_current as rc
import pickle as pkl
import sys


def dlag_Norms ( lag, pers_type, b_num, normalize='',rmv='',norm=2 ):
    cell_Abbrs = ['n4','n5','n11','n13','n14','o5','o9','o10','o12']
    normMap = {}
    for cAbbr in cell_Abbrs:
        print cAbbr
        dvec = dlag_vec ( cAbbr, lag, pers_type, b_num, 
                                            0, 2, normalize, 'bin',rmv)
        normMap[cAbbr] = get_Norm (dvec, norm)
    return normMap
    
def dlag_Norms_More ( lag, pers_type, b_num, normalize='',rmv='',norm=2 ):
    cell_Abbrs = ['n3','n4','n6','n9','o7','o8','o9']
    #cell_Abbrs = ['n3','n4','n9','o7','o8','o9']
    normMap = {}
    for cAbbr in cell_Abbrs:
        print cAbbr
        dvec = dlag_vec_Block ( cAbbr, lag, pers_type, b_num, 
                                            0, 2, normalize, 'bin',rmv)
        normMap[cAbbr] = get_Norm (dvec, norm)
    return normMap
    
def dlag_Norms_All ( lag, pers_type, b_num, lb=0,ub=2,normalize='',rmv='',norm=2 ):
    cell_Abbrs = ['n4','n5','n11','n13','n14','o5','o9','o10','o12']
    cell_Abbrs_More = ['n3','n4','n6','n9','o7','o8','o9']
    normMap = {}
    for cAbbr in cell_Abbrs:
        print cAbbr
        dvec = dlag_vec ( cAbbr, lag, pers_type, b_num, 
                                            lb, ub, normalize, 'bin',rmv)
        normMap[cAbbr] = get_Norm (dvec, norm)
    for cAbbr in cell_Abbrs_More:
        print cAbbr
        dvec = dlag_vec_Block ( cAbbr, lag, pers_type, b_num, 
                                            lb, ub, normalize, 'bin',rmv)
        if cAbbr = 'o9':
            normMap['old_9'] = get_Norm (dvec, norm)
        else:
            normMap[cAbbr] = get_Norm (dvec, norm)

    return normMap

def dlag_Avg ( lag, pers_type, b_num, normalize='',rmv='',norm=2 ):
    cell_Abbrs = ['n4','n5','n11','n13','n14','o5','o9','o10','o12']
    avgMap = {}
    for cAbbr in cell_Abbrs:
        print cAbbr
        dvec = dlag_vec ( cAbbr, lag, pers_type, b_num, 
                         0, 1, normalize, 'bin',rmv,'/data/tOutAvg.txt')
        thisMap = {}
        thisMap['mean'] = numpy.mean(dvec)
        thisMap['median'] = numpy.median(dvec)
        avgMap[cAbbr] = thisMap
    return avgMap

def lagH_Norm ( lag, b_num, normalize='',rmv='',norm=2 ):
    cell_Abbrs = ['n4','n5','n11','n13','n14','o5','o9','o10','o12']
    normMap = {}
    for cAbbr in cell_Abbrs:
        print cAbbr
        dvec = lag_h ( cAbbr, lag, b_num, normalize, rmv, 'bin',2,0,1)
        normMap[cAbbr] = get_Norm (dvec, norm)
    return normMap

def dlag_vec (cAbbr, lag, type, b_num,lb,ub,normalize='', out_type='',rmv='',
              tempOut='/data/tOut.txt', save='NO', output='NIL'):
    """
        Compute lag vector
        Arguments are similar as above
        ex cAbbr = 'n11', 'o9', etc
        BUT hanldes the case '/data/PerseusOutput/original/2d_sparse/New/' (i.e. All in one file)        type = 'bd' or 'wass'
        Except for lag, which is the difference variable for computing v_lag ex. lag = 25
        lb serves as lower_bound for sigma and as lower bound on persistence values for binning
        ub serves as upper_bound for sigma and as upper bound on binning (per_bin)
        out_type is outlier type, ex ='sigma' or 'bin'
    """
    cnames = pkl.load(open('persCellDict.pkl','r'))
    fnames = pkl.load(open('fileNames.pkl','r'))
    cell = cnames[cAbbr]
    files = fnames[cAbbr]
    cellFrames = []
    cellDict = {}
    dlist = os.listdir(cell)
    for f in dlist:
        if f.endswith(str(b_num)+'.txt'):    #correct betti num
            if f.split('-')[0] == files.split('/')[-2]:   #correct file type
                cellPath = cell + f
                fname = files + f[:-6] + '.npy'
                cellFrames.append(cellPath)
                cellDict[f.rstrip('.txt')] = fname
    cellFrames.sort(key=natural_key)
    data = numpy.zeros(len(cellFrames) - (lag-1))#allocate length based on lag
    cellStack = []
    genStack = []
    for ind, g in enumerate(cellFrames):
        #g_gens = rc . get_gens_bin ( g, 1, cellDict[g] )
        gstr = g.split('/')
        g_key = gstr[-1].rstrip('.txt')
        g_gens = rc . get_outlier_gens (g, lb, ub, out_type, rmv, cellDict[g_key])
        #create short name of form '/data/$cell_name'
        gshort = '/'+ gstr[1] + '/' + gstr[-1]
        if normalize:
            maxL = get_Max( cellDict[g_key], 1 )
            g_gens = [(float(gen[0])/maxL,float(gen[-1])/maxL) for gen in g_gens]
        if ind >= lag:
            fshort = cellStack.pop(0)
            f_gens = genStack.pop(0)
            temp_files (fshort, gshort, f_gens, g_gens)
            f_key = fshort.split('/')[-1].rstrip('.txt')
            data[ind-lag] = persistence_distance(type,fshort,gshort,
                                                 cellDict[f_key],cellDict[g_key],tempOut)
            delete_files(fshort,gshort)
        cellStack.append(gshort)
        genStack .append(g_gens)
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
    maxl = get_Max(fname1, 1)
    maxl2 = get_Max(fname2, 1)
    #    print '\n...computing ' + type + ' distance..'
    #print persFile1
    #print persFile2
    if type.startswith('b'):
        call(['/home/kellys/Dropbox/WM/KellySpendlove/metric/kel_metric/bottleneck/bd', persFile1, persFile2, str(maxl), str(maxl2), tempOut])
    else:
        call(['/home/kellys/Dropbox/WM/KellySpendlove/metric/kel_metric/wasserstein/src/wass', persFile1, persFile2, str(p), str(maxl), str(maxl2), tempOut])
    with open(tempOut, 'r') as f:
        d = float(f.read())
    os.remove(tempOut)
    return d

def temp_files (fname, gname, fgens, ggens):
    with open(fname, 'w') as fh:
        rows = map(lambda g: str(g[0])+' '+str(g[-1])+'\n', fgens)
        fh.writelines(rows)
    with open(gname,'w') as fh2:
        rows2 = map(lambda g: str(g[0])+' '+str(g[-1])+'\n', ggens)
        fh2.writelines(rows2)

def delete_files (fname1, fname2):
    os.remove(fname1)
    os.remove(fname2)

def lag_h (cAbbr, lag, b_num,normalize='',rmv='',out_type='',
           p=2,lb=1,ub=3, save='NO', output='NIL'):
    """
        Compute lag vector
        Arguments are similar as above
        ex cAbbr = 'n11', 'o9', etc
        BUT hanldes the case '/data/PerseusOutput/original/2d_sparse/New/' (i.e. All in one file)        type = 'bd' or 'wass'
        Except for lag, which is the difference variable for computing v_lag ex. lag = 25
        """
    cnames = pkl.load(open('persCellDict.pkl','r'))
    fnames = pkl.load(open('fileNames.pkl','r'))
    cell = cnames[cAbbr]
    files = fnames[cAbbr]
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
    genStack = []
    for ind, g in enumerate(cellFrames):
        g_gens = rc . get_outlier_gens (g, lb, ub, out_type, rmv, cellDict[g])
        if normalize:
            maxL = get_Max( cellDict[g], 1 )
            g_gens = [(float(gen[0])/maxL,float(gen[-1])/maxL) for gen in g_gens]
        if ind >= lag:
            f = cellStack.pop(0)
            f_gens = genStack.pop(0)
            data[ind-lag] = hausdorff_distance ( g_gens, f_gens, p)
        cellStack.append(g)
        genStack . append(g_gens)
    if save=='YES':
        saveMatrix(data,output)
    return data

def hausdorff_distance ( goodGensX, goodGensY, p=2 ):
    #intialize sup
    supX = 0
    for genX in goodGensX:
        #convert generator to array
        x = numpy.array(genX)
        #initialize inf
        infY = sys.maxint
        for genY in goodGensY:
            y = numpy.array(genY)
            #get distance
            temp = get_Norm(x-y,p)
            #update inf
            if temp < infY:
                infY = temp
        #update sup
        if infY > supX:
            supX = infY
    #reverse
    supY = 0
    for genY in goodGensY:
        y = numpy.array(genY)
        infX = sys.maxint
        for genX in goodGensX:
            x = numpy.array(genX)
            temp = get_Norm(x-y,p)
            if temp < infX:
                infX = temp
        if infX > supY:
            supY = infX
    return max( supX , supY )

def dlag_gen_stats ( lag, pers_type, b_num, normalize='',rmv='',norm=2 ):
    cell_Abbrs = ['n4','n5','n11','n13','n14','o5','o9','o10','o12']
    genMap = {}
    for cAbbr in cell_Abbrs:
        print cAbbr
        genMap[cAbbr] = good_gens ( cAbbr, lag, pers_type, b_num, 
                         0, 1, normalize, 'bin',rmv)
    return genMap

def good_gens (cAbbr, lag, type, b_num,lb,ub,normalize='', out_type='',rmv='',
              tempOut='/data/tOut.txt', save='NO', output='NIL'):
    """
        Compute lag vector
        Arguments are similar as above
        ex cAbbr = 'n11', 'o9', etc
        BUT hanldes the case '/data/PerseusOutput/original/2d_sparse/New/' (i.e. All in one file)        type = 'bd' or 'wass'
        Except for lag, which is the difference variable for computing v_lag ex. lag = 25
        lb serves as lower_bound for sigma and as lower bound on persistence values for binning
        ub serves as upper_bound for sigma and as upper bound on binning (per_bin)
        out_type is outlier type, ex ='sigma' or 'bin'
        """
    cnames = pkl.load(open('persCellDict.pkl','r'))
    fnames = pkl.load(open('fileNames.pkl','r'))
    cell = cnames[cAbbr]
    files = fnames[cAbbr]
    cellFrames = []
    cellDict = {}
    dlist = os.listdir(cell)
    for f in dlist:
        if f.endswith(str(b_num)+'.txt'):    #correct betti num
            if f.split('-')[0] == files.split('/')[-2]:   #correct file type
                cellPath = cell + f
                fname = files + f[:-6] + '.npy'
                cellFrames.append(cellPath)
                cellDict[f.rstrip('.txt')] = fname
    cellFrames.sort(key=natural_key)
    numGens = []
    birthAvg = []                
    for ind, g in enumerate(cellFrames):
        gstr = g.split('/')
        g_key = gstr[-1].rstrip('.txt')
        g_gens = rc . get_outlier_gens (g, lb, ub, out_type, rmv, cellDict[g_key])
        #average number of generators
        numGens.append( len(g_gens) )
        #birth times
        births = [g_gens[i][0] for i in xrange(len(g_gens))]
        birthsArr = numpy.array(births)
        birthAvg.append( numpy.mean(birthsArr) )
        #birthAvg.append( numpy.median(birthsArr) ) #use median for handling outliers?
        if normalize:
            maxL = get_Max( cellDict[g_key], 1 )
            g_gens = [(float(gen[0])/maxL,float(gen[-1])/maxL) for gen in g_gens]
            births = [g_gens[i][0] for i in xrange(len(g_gens))]
            birthsArr = numpy.array(births)
            birthAvg = []
            birthAvg.append( numpy.mean(birthsArr) )
            #birthAvg.append( numpy.median(birthsArr) ) #use median for handling outliers?
    numGensArr = numpy.array(numGens)
    birthAvgArr = numpy.array(birthAvg)
    numTuple = ( numpy.mean(numGensArr),numpy.median(numGensArr) )
    birthTuple = (numpy.mean(birthAvg),numpy.median(birthAvg) )
    return numTuple, birthTuple
    
def get_Norm ( data, cast=2 ):
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
    
def dlag_vec_Block (cAbbr, lag, type, b_num,lb,ub,normalize='', out_type='',rmv='',
              tempOut='/data/tOut.txt', save='NO', output='NIL'):
    """
        Compute lag vector
        Arguments are similar as above
        ex cAbbr = 'n11', 'o9', etc
        Except for lag, which is the difference variable for computing v_lag ex. lag = 25
        lb serves as lower_bound for sigma and as lower bound on persistence values for binning
        ub serves as upper_bound for sigma and as upper bound on binning (per_bin)
        out_type is outlier type, ex ='sigma' or 'bin'
    """
    fnames = pkl.load(open('moreCellNames.pkl','r'))#location of .npy
    cnames = pkl.load(open('morePersFiles.pkl','r'))#location of perseus input
    cell = cnames[cAbbr]
    file_npy = fnames[cAbbr]
    cellFrames = []
    cellDict = {}
    dlist = os.listdir(cell)
    npy_matrix = numpy.load(file_npy)
    for f in dlist:
        if f.endswith(str(b_num)+'.txt'):    #correct betti num
            cellPath = cell + f
            #fname = files + f[:-6] + '.npy'
            cellFrames.append(cellPath)
            max = npy_matrix[:,:,int(f.split('_')[-2])].max()+1
            cellDict[f.rstrip('.txt')] = max#max for each cell of file
    cellFrames.sort(key=natural_key)
    data = numpy.zeros(len(cellFrames) - (lag-1))#allocate length based on lag
    cellStack = []
    genStack = []
    for ind, g in enumerate(cellFrames):
        gstr = g.split('/')
        g_key = gstr[-1].rstrip('.txt')
        g_gens = rc . get_gens_bin_Block (g, lb, ub, rmv, cellDict[g_key])
        #create short name of form '/data/$cell_name'
        gshort = '/'+ gstr[1] + '/' + gstr[-1]
        if normalize:
            #maxL = get_Max_Block( file_npy,cellDict[g_key], 1 )
            maxL = cellDict[g_key]
            g_gens = [(float(gen[0])/maxL,float(gen[-1])/maxL) for gen in g_gens]
        if ind >= lag:
            fshort = cellStack.pop(0)
            f_gens = genStack.pop(0)
            temp_files (fshort, gshort, f_gens, g_gens)
            f_key = fshort.split('/')[-1].rstrip('.txt')
            data[ind-lag] = persistence_distance_Block(type,fshort,gshort,
                                                 cellDict[f_key],cellDict[g_key],tempOut)
            delete_files(fshort,gshort)
        cellStack.append(gshort)
        genStack .append(g_gens)
    if save=='YES':
        saveMatrix(data,output)
    return data
    
def persistence_distance_Block (type, persFile1, persFile2, maxl, maxl2,tempOut, p=2):
    """
        Compute persistence distance ('bd' or 'wass')
        Call modified metric code in kel_metric
        Preemptive get max heights from .npy for infinite persistence
        p for pth Wasserstein distance (default of 2)
        fname1, fname2 are indices for 3d npy block
    """
    #maxl = get_Max_Block(npy_file, fname1, 1)
    #maxl2 = get_Max_Block(npy_file, fname1, 1)
    #    print '\n...computing ' + type + ' distance..'
    #print persFile1
    #print persFile2
    if type.startswith('b'):
        call(['/home/kellys/Dropbox/WM/KellySpendlove/metric/kel_metric/bottleneck/bd', persFile1, persFile2, str(maxl), str(maxl2), tempOut])
    else:
        call(['/home/kellys/Dropbox/WM/KellySpendlove/metric/kel_metric/wasserstein/src/wass', persFile1, persFile2, str(p), str(maxl), str(maxl2), tempOut])
    with open(tempOut, 'r') as f:
        d = float(f.read())
    os.remove(tempOut)
    return d
    
def get_Max_Block ( fname, ind, add=0 ):
    """
        Send in add=k to return max height value + k
        For example, setting add=1 works when H_1 
    """
    return numpy.load(fname)[:,:,ind].max()+add

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
"""
    
Module for assorted current functions
    created by kel 6/23/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""

import numpy
import matplotlib.pyplot as plt
import re
import os
import time
import cPickle as pkl
from collections import defaultdict

def get_gens_between (file, epsilon1, epsilon2):
    """
        Epsilon1 is lower bound, epsilon2 is upper bound, i.e.
        epsilon1 < (death - birth) < epsilon2
    """
    with open(file, 'r') as f:
        s = f.read()
    f.close()
    goodGens = []
    #split up generators
    stringGens = s.split('\n')
    stringGens.remove('')
    gens = []
    #parse generators
    for sgen in stringGens:
        gens.append(map(int,sgen.split(' ')))
    for (birth,death) in gens:
        if (death - birth) > epsilon1 and (death-birth) < epsilon2:
            goodGens.append((birth,death))
    return goodGens

def get_gens_between_normed( fname, eps1, eps2, means=False ):
    """
    """
    gens = numpy.loadtxt( fname )
    y1, y2 = normalize_mid_lifespan( gens, eps1, eps2 )
        
    # normalize the generator stack and the midrange band
    gens = normalize( gens )
    goodGens = []
    # now find the normalized midrange generators birth and death time
    for (birth,death) in gens:
        if (death - birth) > y1 and (death-birth) < y2:
            goodGens.append((birth,death))
    if means:
        if goodGens:
            return numpy.asarray( goodGens )
        else:
            return None
    else:
        return goodGens

def normalize_mid_lifespan( gens, eps0, eps1 ):
    """
    Normalize the midrange band. Basically, take

    f: (eps0,eps1) --> [0,1].

    Returns f(eps0,eps1) = (y0,y1)
    """
    eps0 = float( eps0 )
    if eps1 == -1:
        eps1 = float( gens.max() )
    else:
        eps1 = float( eps1 )
    delta = gens.max()
    return eps0 / delta, eps1 / delta


def normalize(arr, imin=0, imax=1, dmin=None, dmax=None):
    """
    Normalize 'arr', in-place. (Stolen from stack
    overload. Surprised numpy doesn't have a built-in normalize
    function.)

    (imin, max) -- desired range of normalization

    dmin and dmax -- used if the array does not include all of the
    values. For example, birth time may not include the minimum and
    maximum values. In this case, 0 and max_height are passed to the
    function.
    """
    arr = arr.astype( float )
    if dmin is None:
        dmin = 0 #arr.min()
    if dmax is None:
        dmax = arr.max()
    arr -= dmin
    arr *= (imax - imin)
    arr /= (dmax - dmin)
    arr += imin
    return arr
    
def get_mean_sigma ( fname ):
    """
        Function to get mean and standard deviation
    """
    tsArr = get_ts ( fname )
    mean = tsArr.mean()
    std = tsArr.std()
    return (mean, std)

def get_ts ( file, data='', lb=0, ub=None ):
    """
        Get time series of persistence values for file
        (persistence value is the death-birth integer for each generator)
        file is Perseus output file, ex:
        '/data/PerseusData/PerseusOutput/original/2d_sparse/New/new_110125/
        new_110125-concatenated-ASCII_1000_1.txt'
        i.e.
        '.._xxxx_$betti_num_1.txt'

    ** No error checking that ub argument makes sense **
    """
    gens = get_gens ( file, data )
    if not ub:
        # just set the upper bound to infinity
        ub = numpy.infty
    ts = [(death-birth) for (birth,death) in gens
          if (death-birth) > lb and death < ub]
    tsArr = numpy.array(ts)
    return tsArr

def get_midrange_ts( fdir, lb, betti=1, sname=None, mean=False, infinite=False ):
    """
    fdir -- directory containing persistence diagrams

    lb -- lower bound for the midrange window

    betti -- 0,1,2: which betti number to extract generators for

    sname -- full path for directory to save file
    """
    dlist = os.listdir( fdir )
    diag_list = []
    # find all diagram files in <fdir> for <betti>
    suffix = '_' + str( betti ) + '.txt'
    for f in dlist:
        if f.endswith( suffix ):
            diag_list.append( f )
            fname = fdir + f
    diag_list.sort( key=natural_key )
    midrange_gens = []
    num_gens = []
    if infinite:
        for d in diag_list:
            # account for the one infinite generator
            num_gens.append( len(get_ts( fdir + d, lb=lb ))-1 )
            midrange_gens.append( get_ts( fdir + d, lb=lb )[:-1] )
    # just truncate the noisy generators, leaving the 'infinite' ones in the list
    else:
        for d in diag_list:
            num_gens.append( len(get_ts( fdir + d, lb=lb )) )
            midrange_gens.append( get_ts( fdir + d, lb=lb ) )
    genarr = midrange_gens #numpy.array( midrange_gens, dtype=numpy.int )
    if mean:
        genarr = numpy.asarray( midrange_gens, dtype=numpy.int )
        return genarr.mean(), genarr.std()
    if sname:
        numpy.savetxt( sname, genarr )
    else:
        return genarr

def get_gens ( file, rmv='',data = ''):
    """
        Get generators function
    """
    with open(file, 'r') as f:
        s = f.read()
    f.close()
    #split up generators
    stringGens = s.split('\n')
    stringGens.remove('')
    gens = [map(int,sgen.split(' ')) for sgen in stringGens]
    if data:
        maxHeight = get_Max( data, 1)
        for gen in gens:
            if gen[-1] == -1:
                gen[-1] == maxHeight
    if rmv:
        remove_inf( gens )
    return gens

def get_gens_bin_Block (file, lb=0, per_bin=3,rmv='',max = ''):
    """
        Returns outlying generators whose persistence values
        are unique (=1) (or beneath per_bin requirement)
    """
    gens = get_gens_Block (file,rmv,max )
    pkMap = {}
    pGenMap = {}
    goodGens = []
    for (birth,death) in gens:
        pers = death-birth
        if pers in pkMap:
            pkMap[pers] = pkMap[pers]+1
            pGenMap[pers].append( (birth,death) )
        else:
            pkMap[pers] = 1
            pGenMap[pers] = [(birth, death)]
    #get the lowest persistence value such that all larger 
    #pers values have k <= 2
    lb_pers = 0
    for p in sorted(pkMap.iterkeys()):
        if pkMap[p] > per_bin:
            lb_pers = p
    for p, k in pkMap.iteritems():
        #if k <= 2 and p >= lb_pers: #madalena suggestion
        if p >= lb_pers:
            goodGens.extend (pGenMap[p])
    #print len(goodGens)
    return goodGens

def get_gens_Block ( file, ind, rmv='',max = ''):
    """
        Get generators function
    """
    with open(file, 'r') as f:
        s = f.read()
    f.close()
    #split up generators
    stringGens = s.split('\n')
    stringGens.remove('')
    gens = [map(int,sgen.split(' ')) for sgen in stringGens]
    #print gens
    if max:
        maxHeight = max#get_Max_Block( data, ind, 1)
        for gen in gens:
            if gen[-1] == -1:
                gen[-1] == maxHeight
    if rmv:
        remove_inf( gens )
    return gens
    
def remove_inf_Block ( gens ):
    max = 0
    maxInd = []
    for ind in xrange(len(gens)):
        pers = gens[ind][1] - gens[ind][0]
        if pers > max:
            maxInd = []
            maxInd.append(ind)
        if pers == max:
            #add ind to remove minus number to be removed ahead of it
            maxInd.append( ind-len(maxInd) )
        max = pers #update the current maximum
    for ind in maxInd:
        gens.pop( ind )

def remove_inf ( gens ):
    max = 0
    maxInd = []
    for ind in xrange(len(gens)):
        pers = gens[ind][1] - gens[ind][0]
        if pers > max:
            maxInd = []
            maxInd.append(ind)
        if pers == max:
            #add ind to remove minus number to be removed ahead of it
            maxInd.append( ind-len(maxInd) )
        max = pers
    for ind in maxInd:
        gens.pop( ind )

def get_gens_normalize ( file, data = ''):
    """
        Get generators function
        
        this function is currently not used as techniques for
        retrieving outlying generators is independent
        of normalization (so far) - kel, 7/2/2012
        """
    with open(file, 'r') as f:
        s = f.read()
    f.close()
    #split up generators
    stringGens = s.split('\n')
    stringGens.remove('')
    maxHeight = get_Max( data, 1)
    gens = [map(int,sgen.split(' ')) for sgen in stringGens]
    for gen in gens:
        if gen[-1] == -1:
            gen[-1] == maxHeight
    n_gens = [(float(gen[0])/maxHeight,float(gen[-1])/maxHeight) for gen in gens]
    return n_gens

def get_gens_bin (file, lb=0, per_bin=3,rmv='',data = ''):
    """
        Returns outlying generators whose persistence values
        are unique (=1) (or beneath per_bin requirement)
    """
    gens = get_gens (file,rmv,data )
    pkMap = {}
    pGenMap = {}
    goodGens = []
    for (birth,death) in gens:
        pers = death-birth
        if pers in pkMap:
            pkMap[pers] = pkMap[pers]+1
            pGenMap[pers].append( (birth,death) )
        else:
            pkMap[pers] = 1
            pGenMap[pers] = [(birth, death)]
    #get the lowest persistence value such that all larger 
    #pers values have k <= 2
    lb_pers = 0
    for p in sorted(pkMap.iterkeys()):
        if pkMap[p] > per_bin:
            lb_pers = p
    for p, k in pkMap.iteritems():
        #if k <= 2 and p >= lb_pers: #madalena suggestion
        if p >= lb_pers:
            goodGens.extend (pGenMap[p])
    return goodGens

def get_gens_sigma (file, lb=1.5, ub=6,rmv='',data = ''):
    """
        Epsilon1 is lower bound, epsilon2 is upper bound, i.e.
        epsilon1 < (death - birth) < epsilon2
        """
    gens = get_gens(file, rmv,data)
    ts = [(death-birth) for (birth,death) in gens]
    tsArr = numpy.array(ts)
    mean,std = tsArr.mean(), tsArr.std()
    goodGens = []
    for (birth, death) in gens:
            check = abs((death-birth)-mean)/std
            if check > lb and check < ub:
                goodGens.append((birth,death))
    return goodGens

def get_outlier_gens (file, lb, ub, out_type='bin', rmv='',data=''):
    """
        - lb serves as overload variable for both
        lower bound on sigma and persistence value 
        for sigma, bin methods, respectively
        - type serves as flag for function desired
        - data variable is only required if user wants to 
        replace the -1 inf pers values by max+1
        - rm is a flag to remove persistence values from time series
    """
    if out_type.startswith('b'):
        goodGens = get_gens_bin (file,lb,ub,rmv,data)
    elif out_type.startswith('w'): # window
        goodGens = get_gens_between( file, lb, ub )
    else:
        goodGens = get_gens_sigma (file,lb,ub,rmv,data)
    if len(goodGens) == 0:
        print 'FLAG ------- 0 Generators!!!'
        print file
    return goodGens

def get_Max ( fname, add=1 ):
    """
        Send in add=k to return max height value + k
        For example, setting add=1 works when H_1 
        """
    try:
        return numpy.load(fname).max()+add
    except IOError:
        return int( numpy.loadtxt(fname).max() + add )

def get_Max_Block ( fname, ind, add=1 ):
    """
        Send in add=k to return max height value + k
        For example, setting add=1 works when H_1 
        """
    return numpy.load(fname)[:,:,ind].max()+add

def get_gens_folder ( fdir, betti_num, epsilon1=0, epsilon2=0):
    """
        get number of good generators
        Lower bound epsilon1 upper bound epsilon2, i.e.
        eps1 < (death - birth) < eps2
    """
    if not fdir.endswith('/'):
        fdir+='/'
    if os.path.isdir(fdir):
        dlist = os.listdir(fdir)
        files = []
        for f in dlist:
            if f.endswith(str(betti_num)+'.txt'):
                files.append(fdir+f)
    else:
        print 'Error: input is not a directory'
    files.sort(key=natural_key)
    cell_gens = []
    for file in files:
        #CURRENTLY ONLY NUMBER OF GENERATORS
        cell_gens.append(len(get_gens_sigma(file,epsilon1)))
        #print file     
    #RETURN AVERAGE NUMBER ACROSS FRAMES
    return float(sum(cell_gens))/float(len(cell_gens)), cell_gens

def get_gens_cell (cellFolder, betti_num, eps1=0, eps2=0):
    """
        get all generators for a cell
    """
    cat = cellFolder.split('/')[-1].lower()
    if not cellFolder.endswith('/'):
        cellFolder+='/'
    if os.path.isdir(cellFolder):
        dlist = os.listdir(cellFolder)
        cells = []
        for f in dlist:
            if os.path.isdir(cellFolder+f) and f.startswith(cat):
                cells.append(cellFolder+f)
    cells.sort(key=natural_key)
    all_gens = []
    #ALL GENS IS CURRENTLY LIST OF AVERAGES FOR GIVEN CELL
    for cell in cells:
        all_gens.append((cell.split('/')[-1],get_gens_folder(cell,betti_num,eps1,eps2)))
    return all_gens

def compute_generator_stats( prefix, cell_type='new', cell_list=None, betti=1, lb=1, ub=10 ):
    """
    """
    if not cell_list:
        new_list = [ 'new_110125', 'new_130125', 'new_140125', 'new_40125', 'new_50125']
        old_list = [ 'old_100125', 'old_120125', 'old_50125', 'old_90125' ]
    if cell_type == 'new':
        prefix += 'New/'
        cell_list = new_list
    elif cell_type == 'old':
        prefix += 'Old/'
        cell_list = old_list
    else:
        print "Unknown cell type!"
        exit(0)
    # d['cell'] = dict{ k : ( mean, std ) }
    outlier_stats = {}
    # loops over all cells, then record the mean and std for each cell
    for cell in cell_list:
        print "Computing outlier stats for cell", cell
        outlier_stats[ cell ] = {}
        for k in numpy.arange( 1, 5.5, 0.5 ): # excludes endpoint
            outlier_mean, outliers = get_gens_folder( prefix+cell, betti,
                                                      epsilon1=k*lb, epsilon2=ub )
            outlier_stats[ cell ][k] = ( outlier_mean,
                                         numpy.std( outliers ) )
    with open( 'outlier_stats_'+cell_type+'.pkl', 'w' ) as fh:
        pkl.dump( outlier_stats, fh, protocol=-1 )
    return outlier_stats

def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
    
    
def get_short_range_gens (file,per_bin=3,rmv='Y',data = ''):
    """
        Returns outlying generators whose persistence values
        are unique (=1) (or beneath per_bin requirement)
    """
    gens = get_gens (file,rmv,data )
    pkMap = {}
    pGenMap = {}
    shortGens = []
    for (birth,death) in gens:
        pers = death-birth
        if pers in pkMap:
            pkMap[pers] = pkMap[pers]+1
            pGenMap[pers].append( (birth,death) )
        else:
            pkMap[pers] = 1
            pGenMap[pers] = [(birth, death)]
    #get the lowest persistence value such that all larger 
    #pers values have k <= 2
    lb_pers = 0
    for p in sorted(pkMap.iterkeys()):
        if pkMap[p] > per_bin:
            lb_pers = p
    for p, k in pkMap.iteritems():
        #if k <= 2 and p >= lb_pers: #madalena suggestion
        if p < lb_pers:
            shortGens.extend (pGenMap[p])
    return shortGens

def plot_sgen ( b_num ):
    cell_Abbrs = ['n4','n5','n11','n13','n14','o5','o9','o10','o12']
    cnames = pkl.load(open('persCellDict.pkl','r'))
    fnames = pkl.load(open('fileNames.pkl','r'))
    cellMap = {}
    for cAbbr in cell_Abbrs:
        dlist = os.listdir ( fnames[cAbbr] )
        frames = []
        gr_Num = 0
        for f in dlist:
            if f.endswith('.npy'):
                frames.append ( f.rstrip('.npy')+'_'+str(b_num)+'.txt' )
                if not gr_Num:
                    data = numpy.load(fnames[cAbbr]+f )
                    gr_data = data[data>0]
                    gr_Num = len(gr_data)
        lenSGen = []
        for frame in frames:
            shortGens = get_short_range_gens ( cnames[cAbbr]+frame )
            lenSGen . append ( len(shortGens) )
        lsgen = numpy.array(lenSGen)
        lsNum = lsgen . mean()
        cellMap[cAbbr] = [gr_Num, lsNum]
    x = []
    y = []
    for c,numList in cellMap.iteritems():
        x . append (numList[0])
        y . append (numList[-1])
    plt.scatter(x,y)
    plt.show()

def dir_list( fdir, betti=1 ):
    """
    Returns a list of Perseus output files for given betti #.
    """
    dlist = os.listdir( fdir )
    theFiles = [ fdir+f for f in dlist if f.endswith( '_'+str(betti)+'.txt' ) ]
    theFiles.sort( key=natural_key )
    return theFiles


if __name__ == "__main__":

    #prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
    prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
    newlist = ['new_10', 'new_110125', 'new_130125', 'new_140125', 'new_3',
               'new_4', 'new_40125', 'new_50125', 'new_6', 'new_60125', 'new_9']
    oldlist = ['old_100125', 'old_120125', 'old_15', 'old_2', 'old_4000', 'old_4001',
               'old_5',  'old_50125',  'old_6',  'old_7',  'old_8',  'old_9',  'old_90125']
    #cells = [ prefix + c + '/' for c in newlist ]
    cells = [ prefix + c + '/' for c in oldlist ]
    frames = [ dir_list( c ) for c in cells ]
    # fig, ts = plot_hist_stack( frames, left_xlim=0.2, right_xlim=0.6, normed=False,
    #                            cutoff=0.2, ts_max=1000, skip=20, log=True )

    # lower bound
    lb = [55, 60]
    old_ts = {}
    for x in lb:
        for cell_dir in cells:
            print "Getting midrange gens for ", cell_dir
            old_ts[ cell_dir ] = get_midrange_ts( cell_dir, x )
        with open( './timeseries/old_gen_ts_lb'+str(x)+'.pkl', 'w' ) as fh:
            pkl.dump( old_ts, fh )

    # lb = [60, 73]
    # new_ts = {}
    # for x in lb:
    #     for cell_dir in cells:
    #         print "Getting midrange gens for ", cell_dir
    #         new_ts[ cell_dir ] = get_midrange_ts( cell_dir, x )
    #     with open( './timeseries/new_gen_ts_lb'+str(x)+'.pkl', 'w' ) as fh:
    #         pkl.dump( new_ts, fh )

    # plot_hist_cut_axis( frames, normed=False, log=True,
    #                     left_xlim=0.2, right_xlim=0.6,
    #                     cutoff=0.2, ts_max=1000, skip=20 )

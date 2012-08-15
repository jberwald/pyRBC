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
    [gens.append(map(int,sgen.split(' '))) for sgen in stringGens]
    for (birth,death) in gens:
        if (death - birth) > epsilon1 and (death-birth) < epsilon2:
            goodGens.append((birth,death))
    return goodGens

def get_mean_sigma ( file ):
    """
        Function to get mean and standard deviation
    """
    tsArr = get_ts ( file )
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
    #return ts #return list
    return tsArr

def get_midrange_ts( fdir, lb, betti=1, sname=None ):
    """
    fdir -- directory containing persistence diagrams

    lb -- lower bound for the midrange window

    betti -- 0,1,2: which betti number to extract generators for

    sname -- full path for directory to save file
    """
    dlist = os.listdir( fdir )
    diag_list = []
    # find all diagram files in <fdir> for <betti>
    for f in dlist:
        if f.endswith( '_'+str( betti )+'.txt'):
            diag_list.append( f )
            fname = fdir + f
    diag_list.sort( key=natural_key )
    midrange_gens = []
    num_gens = []
    for d in diag_list:
        # account for the one infinite generator
        num_gens.append( len(get_ts( fdir + d, lb=lb ))-1 )
        midrange_gens.append( get_ts( fdir + d, lb=lb )[:-1] )
    genarr = numpy.array( midrange_gens, dtype=numpy.int )
    ng
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
    return numpy.load(fname).max()+add

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
    
def plot_hist( fname, color='blue', normed=False, fontsize=20 ):
    """
    Plot a histogram of generator lifespan along the diagonal.

    fname -- full path to perseus output file.
    """
    ts = get_ts ( fname )
    # the (almost) infinite generator overwhelms the plot
    ts = ts[:-1]
    fig = plt.figure()
    ax = fig.gca()
    # the histogram of the data
    n, bins, patches = ax.hist( ts, bins=ts.max()-ts.min(),
                                normed=normed, facecolor=color, alpha=0.75)
    xticks = [ int( tk ) for tk in ax.get_xticks() ]
    yticks = [ int( tk ) for tk in ax.get_yticks() ]
    ax.set_xticklabels( xticks, fontsize=fontsize )
    ax.set_yticklabels( yticks, fontsize=fontsize )

    #ax.set_title( r'Distribution of generator lifespans along diagonal' )
    ax.set_xlabel( r"Lifespan (death-birth)", fontsize=fontsize )
    ax.set_ylabel( r"Number of generators ($\beta_1$)", fontsize=fontsize )
    ax.grid( True )
    plt.show()
    return fig

def dir_list( fdir, betti=1 ):
    """
    Returns a list of Perseus output files for given betti #.
    """
    dlist = os.listdir( fdir )
    theFiles = [ fdir+f for f in dlist if f.endswith( '_'+str(betti)+'.txt' ) ]
    theFiles.sort( key=natural_key )
    return theFiles

def plot_hist_colors( cell, color='blue',
                      normed=False, fontsize=20,
                      threshold=50, cell_type='New',
                      show_plot=False, log=True):
    """
    Plot a histogram of generator lifespan along the diagonal. This
    allows more control over bin color.

    cell -- full path to perseus output file.
    """
    from matplotlib.ticker import ScalarFormatter

    fig = plt.figure( dpi=80)
    ax = fig.gca()

    # cell is a list of frames
    if hasattr( cell, 'pop' ):
        # create an array with first entry in cell, then extend in
        # loop
        ts = get_ts ( cell[0] )
        #ts = ts[:-1]
        for f in cell[1:]:
            new_ts = get_ts ( f )
            # the (almost) infinite generator overwhelms the plot
            #new_ts = new_ts[:-1]
            # extend last ts of generators by new_ts
            ts = numpy.hstack( ( ts, new_ts ) )
            # the histogram of the data
            #vals = vals[1:]
        # create two identical histograms 
        n, bins, patches = ax.hist( ts, bins=(ts.max()-ts.min()),\
                                    normed=normed, facecolor=color,\
                                    alpha=0.75, histtype='bar', log=log )
    # or cell is just a single frame
    else:
        ts = get_ts ( cell )
        # the (almost) infinite generator overwhelms the plot
        ts = ts[:-1]
        # the histogram of the data
        n, bins, patches = ax.hist( ts, bins=ts.max()-ts.min(),\
                                    normed=normed, facecolor=color,\
                                    alpha=0.75 )#  histtype='step')

# # we need to normalize the data to 0..1 for the full
# # range of the colormap
# fracs = N.astype(float)/N.max()
# norm = colors.normalize(fracs.min(), fracs.max())

    theColors = [ 'b', 'r' ]
    for num_gens, patch in zip( bins[:-1], patches ):
        if num_gens > threshold:
            color = theColors[ 1 ]
            # patch.set_width( 1.0 )
        else:
            color = theColors[ 0 ]
            #patch.set_width( 1.0 )
        patch.set_facecolor( color )
    # set some axis attributes
    xticks = [ int( tk ) for tk in ax.get_xticks() ]
    yticks = [ int( tk ) for tk in ax.get_yticks() ]
    ax.set_xticklabels( xticks, fontsize=fontsize )
    ax.set_yticklabels( yticks, fontsize=fontsize )

    #ax.set_title( r''+str( cell_type ) + ' cell', fontsize=24 )
    ax.set_xlabel( r"Lifespan (death-birth)", fontsize=fontsize )
    ax.set_ylabel( r"$\beta_1$", fontsize=fontsize )
    #ax.ticklabel_format( style='sci', axis='y' )  # scientific notation

    sf = ScalarFormatter()
    sf.set_scientific( True )

    ax.grid( True )
    if show_plot:
        plt.show()
    return fig, ts


def plot_midrange_ts( new_file, old_file, skip=10, fontsize=20, lines=None ):
    """
    Plots two times series of midrange generators.
    """
    ts_new = numpy.loadtxt( new_file )
    ts_old = numpy.loadtxt( old_file )

    fig = plt.figure( dpi=160, figsize=([10,4]) )
    ax = fig.gca()
    ax.plot( ts_new[::10], 'b-', linewidth=2 )
    ax.plot( ts_old[::10], 'r-', linewidth=2 )
    # plot vertical lines to point out location of sublevel sets in ts
    if lines:
        vmin, vmax = ax.get_ylim()
        smidge = 0.1
        print vmin, vmax
        for line in lines:
            v = int( line/float(skip) )
            ax.vlines( v, vmin+smidge, vmax, linestyle='dashed', linewidth=2 )


    # set some axis attributes; account for skip when setting tick marks
    xticks = [ skip*int( tk ) for tk in ax.get_xticks() ]
    yticks = [ int( tk ) for tk in ax.get_yticks() ]
    ax.set_xticklabels( xticks, fontsize=fontsize )
    ax.set_yticklabels( yticks, fontsize=fontsize )
    ax.set_xlabel( r"time", fontsize=fontsize )
    ax.set_ylabel( r"# of midrange generators", fontsize=fontsize )
    fig.show()
    return fig

def plot_hist_cut_axis( cell, color='blue',
                        normed=False, fontsize=20,
                        threshold=50, cell_type='New',
                        show_plot=False, log=False,
                        left_xlim=300, right_xlim=1600):
    """
    Plot a histogram of generator lifespan along the diagonal. This
    allows more control over bin color.

    cell -- full path to perseus output file.
    """
    from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes
    from mpl_toolkits.axes_grid.inset_locator import mark_inset
    
    # cell is a list of frames
    # create an array with first entry in cell, then extend in
    # loop
    ts = get_ts ( cell[0] )
    #    ts = ts[:-1]
    for f in cell[1:]:
        new_ts = get_ts ( f )
        # the (almost) infinite generator overwhelms the plot
        #      new_ts = new_ts[:-1]
        # extend last ts of generators by new_ts
        ts = numpy.hstack( ( ts, new_ts ) )
        # the histogram of the data
        #vals = vals[1:]
 

        # This method works with mpl 0.99 (plt.subplot() works with >1.1.0 )
        fig = plt.figure( dpi=160 )
        ax = fig.add_subplot( 121 ) 
        ax2 = fig.add_subplot( 122 )

        #### DEBUGGING ####
        ts = ts[:100]

        # create two identical histograms 
        n, bins, patches = ax.hist( ts, bins=(ts.max()-ts.min()),\
                                    normed=normed, facecolor=color,\
                                    alpha=0.75, histtype='bar', log=log, rwidth=20 )
        n, bins, patches = ax2.hist( ts, bins=(ts.max()-ts.min()),\
                                     normed=normed, facecolor=color,\
                                     alpha=0.75, histtype='bar', log=log, rwidth=20 )

        # make the inset
        # options of loc: BEST, UR, UL, LL, LR, R, CL, CR, LC, UC, C = range(11)
        #axins = zoomed_inset_axes(ax, 10, loc=10 ) # ( axes, zoom power, location )
        axins = inset_axes( ax, width='70%', height=1., loc = 10 )
        axins.hist( ts, bins=(ts.max()-ts.min()),\
                    normed=normed, facecolor=color,\
                    alpha=0.75, histtype='bar', log=log, rwidth=20 )
        # set box position by hand

        # pos = [left, bottom, width, height]
        # axins.set_position

        # hide spines between axes
        ax.spines['right'].set_visible(False)
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        # the tail side
        ax2.spines['left'].set_visible(False)
        ax2.xaxis.tick_bottom()
        ax2.yaxis.tick_right()
        #ax2.set_yticks( [] )
        #ax.tick_params(labeltop='off') # don't put tick labels at the top


        d = .015 # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        # remember, plot takes list of x's, followed by list of y's
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((1-d,1+d),(-d,+d), **kwargs)      # bottom-right diagonal
        ax.plot((1-d,1+d),(1-d,1+d), **kwargs)    # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the right axes
        ax2.plot((-d,+d),(-d,+d), **kwargs)   # bottom-left diagonal
        ax2.plot((-d,+d),(1-d,1+d), **kwargs) # top-left diagonal

       # zoom-in / limit the view to different portions of the data
        ax.set_xlim( 0., left_xlim ) # most of the data
        ax.set_ylim( 0., n.max()+5 )
        ax2.set_xlim( right_xlim, ts.max()+5) # outliers/inf gens
        ax2.set_ylim( 0., 60. )  # zoom in to see the inf gens

        # sub region of the original image
        yscale = n.max() * (1./60)
        xmax = 250.
        x1, x2, y1, y2 = 25, xmax, 0, yscale
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_aspect( xmax/yscale )
        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    fig.subplots_adjust( wspace=0.1 )
    plt.show()
    
    return fig, ts

def plot_hist_stack( cells, color='blue',
                     normed=False, fontsize=20,
                     threshold=50, cell_type='New',
                     show_plot=False, log=False,
                     left_xlim=300, right_xlim=1600,
                     nbins=1000.):
    """
    cells -- list of cells (full paths to) whose generator lifespans
    we want to stack in a single histogram.
    """
    from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes
    from mpl_toolkits.axes_grid.inset_locator import mark_inset
    # ==========  ========
    # character   color
    # ==========  ========
    # 'b'         blue
    # 'g'         green
    # 'r'         red
    # 'c'         cyan
    # 'm'         magenta
    # 'y'         yellow
    # 'k'         black
    # 'w'         white
    # ==========  ========
    colors = ['b','g','r','c','m','y']
    nbins = float( nbins )
    # cell is a list of frames
    # create an array with first entry in cell, then extend in
    # loop

    # This method works with mpl 0.99 (plt.subplot() works with >1.1.0 )
    fig = plt.figure( dpi=160 )
    ax = fig.add_subplot( 121 ) 
    ax2 = fig.add_subplot( 122 )

    # cell == cell directory containing persistence diagrams
    old_ts = None
    for cell in cells:
        ts = get_ts ( cell[0] )
        # concatenate all of the diagram lifespans into one timeseries
        for f in cell[1:]:
            new_ts = get_ts ( f )
            # extend last ts of generators by new_ts
            ts = numpy.hstack( ( ts, new_ts ) )
        # compute the histogram for ts
        h, bins = numpy.histogram( ts, bins=(ts.max()-ts.min())/nbins )
        ax.bar( ts, width=10, bottom=old_ts )
        ax2.bar( ts, width=10, bottom=old_ts )
    
       
    # pos = [left, bottom, width, height]
    # axins.set_position
        
    # hide spines between axes
    ax.spines['right'].set_visible(False)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    # the tail side
    ax2.spines['left'].set_visible(False)
    ax2.xaxis.tick_bottom()
    ax2.yaxis.tick_right()
    #ax2.set_yticks( [] )
    #ax.tick_params(labeltop='off') # don't put tick labels at the top


    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    # remember, plot takes list of x's, followed by list of y's
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d),(-d,+d), **kwargs)      # bottom-right diagonal
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the right axes
    ax2.plot((-d,+d),(-d,+d), **kwargs)   # bottom-left diagonal
    ax2.plot((-d,+d),(1-d,1+d), **kwargs) # top-left diagonal

   # zoom-in / limit the view to different portions of the data
    ax.set_xlim( 0., left_xlim ) # most of the data
    ax.set_ylim( 0., n.max()+5 )
    ax2.set_xlim( right_xlim, ts.max()+5) # outliers/inf gens
    ax2.set_ylim( 0., 60. )  # zoom in to see the inf gens

    # sub region of the original image
    yscale = n.max() * (1./60)
    xmax = 250.
    x1, x2, y1, y2 = 25, xmax, 0, yscale
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # axins.set_aspect( xmax/yscale )
    # # draw a bbox of the region of the inset axes in the parent axes and
    # # connecting lines between the bbox and the inset axes area
    # mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    fig.subplots_adjust( wspace=0.1 )
    plt.show()
    
    return fig, ts

"""
Module for plotting histograms of generator data.
"""
from rbc_current import *
import numpy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
import matplotlib.colors as colors
import timer

slash = '/'

def dict_of_cell_dirs( cell_dirs ):
    """
    cell_dirs -- list of directories (cells) containing cell data.
    """
    # handles last entry in list returned by split
    if cell_dirs[0].endswith( slash ):
        last = -2
    else:
        last = -1
    prefix = slash.join( cell_dirs[0].split( slash )[:last] )
    prefix += slash
    print prefix
    
    cellnames = [ c.split( slash )[last] for c in cell_dirs ]
    print cellnames
    
    cells = dict.fromkeys( cellnames )
    for k in cells.keys():
        cells[k] = dir_list( prefix + k + slash )
    return cells

# def load_birth_times( fname ):
#     with open( fname ) as fh:
#         return pkl.load( fh )

def concat_histograms( ts, width=10 ):
    """
    ts -- dictionary of time series. See concat_timeseries() below.
    """
    all_hist = []
    # set common bin width
    dmax = 0
    for d in ts.itervalues():
        dmax = max
    for i, data in enumerate( ts.itervalues() ):
        thebins = numpy.arange( 0, data.max(), width )
        hist, bins = numpy.histogram( data, bins=thebins )
        all_hist.append( hist )
    return all_hist

def cell_max( fdir, dim=1 ):
    """
    For normalization purposes, find the maximum height over all
    frames for a given cell. 
    """
    if not fdir.endswith( slash ): fdir += slash
    frames = os.listdir( fdir )
    ending = str( dim ) + '.txt'
    frames = [ fdir+f for f in frames if f.endswith( ending ) ]
    the_max = 0
    for frame in frames:
        # from rbc_current
        x = get_Max( frame )
        if x > the_max:
            the_max = x
    return the_max

def all_maxes( cell_list ):
    """
    """
    max_list = [ ( cell, int( cell_max( cell ) ) )
                 for cell in cell_list ]
    return max_list
        

def get_birth_times( cell, eps1, eps2=150, normed=True ):
    """
    cell -- path to directory containing Perseus generator
    file (one for each frame).

    eps1 -- minimum lifespan

    eps2 -- maximum lifespan

    (So [eps1,eps2], stretched along the diagonal, is the band that we
    store generators from.)

    For each frame in cell==cell_dir: find generator lifespans.  Find first
    occurence, \tau, of a midrange generator ( use get_gens_between()
    for this, then peel off the birth time from the first
    (birth,death) pair to get birth time) store \tau.

    NOTE: Depending on eps1, very rarely a list of gens will be
    empty. We treat this as missing data and continue to loop over the
    frames.
    """
    frames = dir_list( cell )
    birth_times = []
    for frame in frames:
        if normed:
            gens = get_gens_between_normed( frame, eps1, eps2 )
        else:
            gens = get_gens_between( frame, eps1, eps2 )
        if not gens:
            continue
        birth_times.append( gens[0][0] )
    return birth_times


def get_mean_midrange( cell, eps1, eps2=150, normed=True):
    """
    cell -- path to directory containing Perseus generator
    file (one for each frame).

    eps1 -- minimum lifespan

    eps2 -- maximum lifespan

    (So [eps1,eps2], stretched along the diagonal, is the band that we
    store generators from.)

    For each frame in cell==cell_dir: find generator lifespans.  Find first
    occurence, \tau, of a midrange generator ( use get_gens_between()
    for this, then peel off the birth time from the first
    (birth,death) pair to get birth time) store \tau.

    NOTE: Depending on eps1, very rarely a list of gens will be
    empty. We treat this as missing data and continue to loop over the
    frames.
    """
    frames = dir_list( cell )
    gen_stats = []
    for frame in frames:
        # get the midrange gen stats for frame (normed or not)
        if normed:
            gstats = get_gens_between_normed( frame, eps1, eps2, means=True )
        else:
            gstats = get_gens_between( frame, eps1, eps2, means=True )
        if not gstats:
            continue
        gen_stats.append( gens )
    return gen_stats

def plot_boxplot( data, vert=1, pa=True, transparent=True ):
    """
    data -- a vector or list of vectors. Can be of various
    lengths. Originally created for vectors consisting of lag 1 norms
    of persistence distances.

    vert -- vert==1 ==> vertical orientation;
            vert==0 ==> horizontal orientation
            (follows boxplot() convention)

    pa -- patch_artist: set to True for solid boxes.

    Produces a boxplot (see pylab doc).
    """
    fig = plt.figure()
    if transparent:
        fig.patch.set_alpha( 0.0 )
    ax = fig.gca()
    ax.boxplot( data, vert=vert, patch_artist=pa )

    if vert == 1:
        ax.set_xlabel( r'Cell type', fontsize=24 )
        ax.set_ylabel( r'Lag 1 norm', fontsize=24 )
        # rename the xlabels
        ax.set_xticklabels( ['Young', 'Old'], fontsize=20 )
        # make the ylabels bigger
        yticks = ax.get_yticks()
        ax.set_yticklabels( [str(int(y)) for y in yticks], fontsize=20 )
    # vert better equal 0
    else:
        ax.set_ylabel( r'Cell type', fontsize=24 )
        ax.set_xlabel( r'Lag 1 norm', fontsize=24 )
        # rename the xlabels
        ax.set_yticklabels( ['Young', 'Old'], fontsize=20 )
        # make the ylabels bigger
        yticks = ax.get_xticks()
        ax.set_xticklabels( [str(int(y)) for y in yticks], fontsize=20 )
    fig.show()
    return fig

def boxplotter( fname='/Users/jberwald/github/local/caja-matematica/pyRBC/data/lag1_all.pkl',
                **kwargs ):
    """
    fname -- Pickled dictionary with lag k distance vector norms.
    """
    with open( '/Users/jberwald/github/local/caja-matematica/pyRBC/data/lag1_all.pkl' ) as fh:
        all_norms = pkl.load( fh )
    old = [ all_norms[k] for k in all_norms if 'o' in k ]
    new = [ all_norms[k] for k in all_norms if 'n' in k ]
    norm_vecs = [ new, old ]
    fig = plot_boxplot( norm_vecs, **kwargs )
    return fig


def plot_hist_birth_times( bt_old=None, bt_new=None,
                           normalize=False, bins=50,
                           transparent=True, **kwargs ):
    """
    Plot histograms of birth times. Provide both to plot both new and
    old on one histogram.

    bt_* -- list of birth times for each cell. 
    """
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection

    fig = plt.figure( figsize=(10,8) )
    if transparent:
        fig.patch.set_alpha( 0.0 )
    ax = fig.gca()
    if bt_old:
        # concatenate data
        old = numpy.hstack( bt_old )
        oldhist, bins, patches = ax.hist( old, bins=bins, hatch='//',
                                          color='r', alpha=0.75, )
        mean_old = old.mean()
        std_old = old.std()
        ax.axvline( mean_old, 
                    color='k', linestyle='dashed', lw=1.5 )
        ax.axvspan( mean_old - std_old,
                    mean_old + std_old, facecolor='r', alpha=0.3)
            # add an arrow -- NOT WORKING YET
        patches = []
        arrow = mpatches.ArrowStyle( "Fancy, head_length=.4, head_width=.4, tail_width=.4") 
        patches.append(arrow)
        # plt.text(pos[0,5], pos[1,5]-0.15, "Arrow", ha="center",
        #          family=font, size=14)
    if bt_new:
        new = numpy.hstack( bt_new )
        newhist, bins, patches = ax.hist( new, bins=bins,
                                          color='b', alpha=0.75 )
        mean_new = new.mean()
        std_new = new.std()
        ax.axvline( mean_new,
                    color='k', linestyle='dashed', lw=1.5 )
        ax.axvspan( mean_new - std_new,
                    mean_new + std_new, facecolor='b', alpha=0.3)

    if not bt_old and not bt_new:
        print "Must provide at least one set of birth times."
        return None
    ax.set_xlabel( r'$\tau$ (normalized)', fontsize=24 )
    ax.set_ylabel( 'Number of frames', fontsize=24 )
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    ax.set_xticklabels( [str( x ) for x in xticks], fontsize=20 )
    ylabels = []
    for y in yticks:
        if y!=0:
            ylabels.append( str( int(y) ) )
        else:
            ylabels.append( ' ' ) # empty 0 on y axis
    ax.set_yticklabels( ylabels, fontsize=20 )
    # ax.set_yticklabels( [str( int(y) ) for y in yticks
    #                      if y!=0 else ], fontsize=20 )

    # This part adds the arrows defined above.
    # collection = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    # ax.add_collection(collection)

    # set axes limits to align with top of vlines for aesthetic
    # reasons
    #ax.set_ylim( 0.0, old_max ) #max( old_max, new_max ) )
    
    
    fig.show()
    return fig, new, old


def plot_hist_all( ts, nbins=50, transparent=True, norm_it=False, **kwargs ):
    """
    ts -- dictionary with values as lifespan timeseries.

    kwargs -- see pylab hist() function

    Returns interpolation functions as well as histogram triple and figure instance.
    """
    from pylab import log1p
    from matplotlib.mlab import stineman_interp

    #data = ts.values()
    data = ts
    fig = plt.figure( figsize=(12,8) )
    if transparent:
        fig.patch.set_alpha( 0.0 )
    # ax = fig.add_subplot( 121 )
    # ax2 = fig.add_subplot( 122 )
    ax2 = fig.add_subplot( 111 )
    # now plot a single cell's histogram on the first axis
    #ax.hist( data, **kwargs )

    xmax = max( [d.max() for d in data] )
    thebins = numpy.linspace(0, xmax, nbins )
    kwargs['bins'] = thebins
    # holds bins counts (n) for each histogram
    all_ny = []
    for d in data:
        n, bins, patches = ax2.hist( d, **kwargs )
        all_ny.append( n )
        # all_bins.append( bins )
    # convert bins counts to a single array to find min
    arr = numpy.array( all_ny ).ravel()
    wy = arr[ numpy.where( arr != 0 )[0] ]
    min_y = min( wy )

    # for plotting average -- these are already log values if
    # log==True in kwargs
    yhist = numpy.array( all_ny, dtype=numpy.float64 )#.ravel()
    avg = yhist.mean( axis=0 )
    err = yhist.std( axis=0 )
    upper = avg + err
    lower = avg - err

    # print "yhist", yhist
    # print "avg", avg
    # print ""

    # print err
    # print ""
    # print lower
    # print upper
    # print ""

    # average value for each histogram bin
    for i, x in enumerate( avg ):
        if x == 0.0:
            avg[i] = 1.0

    # label the axes
    ax2.set_xlabel( 'Generator lifespan', fontsize=20 )
    ax2.set_ylabel( 'Number of generators (log)', fontsize=20 )
    xticks = ax2.get_xticks()
    yticks = ax2.get_yticks()
    ax2.set_xticklabels( [str(int(x)) for x in xticks], fontsize=20 )
    ax2.set_yticklabels( [str(int(y)) for y in yticks], fontsize=20 )

    # now plot the interpolated average distribution here so it is on
    # top of the other stuff
    yp = None
    xi = numpy.linspace( 0, bins[-1],200)
    yi = stineman_interp( xi, bins[:-1], avg, yp )
    # interpolate upper and lower error bars to get envelope
    # COMMENTED OUT BELOW
    upper_yi = stineman_interp( xi, bins[:-1], upper, yp )
    lower_yi = stineman_interp( xi, bins[:-1], lower, yp )
    
    # make sure lower does not go negative since this makes no sense.
    for i,v in enumerate( lower_yi ):
        if v < 1.0:
            lower_yi[i] = 1.0

    # make sure that the plot doesn't get messed up by small values
    # (esp. New cells)
    masked_yi = numpy.ma.masked_less( yi, 1 )

    # plot the interpolation of the avg and the envelope
    ax2.plot( xi, masked_yi, 'r-', lw=3 )
    # ax2.fill_between( xi, lower_yi, upper_yi, #where=masked_yi
    #                   color='r', alpha=0.5, zorder=10 )
    fig.show()

    if norm_it:
        y_max = masked_yi.max()
        print "y_max", y_max
        masked_yi /= y_max

    return xi, masked_yi, lower_yi, upper_yi, fig, (n, bins, patches)

def plot_hist_figure( ts, persfile=None, single=0, norm_it=True, nbins=100 ):
    """
    Plot the histogram figure for the RBC paper.

    ts -- List of arrays of times series of generator lifespans for
    each cell.

    persfile -- full path to persistence file

    single -- Cell to choose from list to compute histogram statistics
    on. Should be the same cell used in <persfile>.

    norm_it -- Toggle whether to return a normalized histogram.
    """
    # compute stats for all cells
    out_all = plot_hist_all( ts, norm_it=norm_it )
    allx = out_all[0]
    ally = out_all[1]

    # compute stats for chosen single cell
    out = plot_hist_all( [ts[single]], norm_it=norm_it )
    nx = out[0]
    ny = out[1]

    # now normalize everything by dividing by total area ( y --> PDF )
    pdf_ally = pdf( allx, ally )
    pdf_ny = pdf( nx, ny )


    print "pdf_ally", ((nx[1:]-nx[:-1]) * pdf_ally[:-1]).sum()
    print "pdf_ny", ((nx[1:]-nx[:-1]) *pdf_ny[:-1]).sum()

    fig = plt.figure()
    ax = fig.gca()
    #    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    #ax.set_aspect( 1 )
 
    ax.plot( nx, pdf_ally, lw=3, c='b', label='Mean, all cells' )
    ax.plot( nx, pdf_ny, lw=3, c='r', label='Mean, single cell' )

    # add a histogram for a single frame
    if persfile:
        ts = numpy.asarray( get_ts( persfile ), dtype=numpy.int )
        n, bins = numpy.histogram( ts, bins=len(nx)-1, range=(nx.min(),nx.max()) )
        ts_pdf = pdf( bins[:-1], n )

        print "ts_pdf", ((bins[1:] - bins[:-1])*ts_pdf).sum()
        width = nx[1]-nx[0]
        ax.bar( bins[:-1], ts_pdf, width=width, color='g', alpha=0.7, label='Single frame distribution' )
        #ax.plot( bins[:-1], ts_pdf, marker='o', ms=6, lw=3, label='Single frame' )
    plt.legend()
    fig.show()
    return fig, bins, ts_pdf  #, n, bins

def pdf( xi, yi ):
    """
    Normalize f(x) = y in terms of probability distribution
    functions. Thus, it should return f such that 

    \int_{min(xi)}^{max(xi)} {f(x_i)*(x_{i+1}-x_{i})   = 1
    """
    dx = xi[1:] - xi[:-1]
    integ = numpy.sum( yi[:-1] * dx )
    thepdf = yi/integ
    return thepdf


def concat_timeseries( cells, ts_max=-1, skip=1, normed=False ):
    """
    cells -- list of cell names (full paths to directories).

    ts_max -- number of frames (so max length of ts). [default=-1, eg. all frames]

    skip -- keep every <skip> frames

    Returns dicitonary, values times series of generator data, keyed
    by cell names
    """
    all_ts = {}
    for c in cells:
        # list of lifespans for each frame in cell
        # concatenate all of the diagrams for cell lifespans into one
        # timeseries
        ts = [ get_ts ( frame )
               for frame in cells[c][:ts_max:skip] ]
        ts = numpy.hstack( ts )
        # convert and normalize, then append to the list of time series
        if normed:
            ts = numpy.asarray( ts, dtype=numpy.float )
            ts /= ts.max()
        else:
            ts = numpy.asarray( ts, dtype=numpy.int )
        all_ts[c]= ts
    return all_ts

def plot_scatter( ts, log=False, cutoff=0.2 ):
    """
    ts -- dictionary of time series data. keys are cell names and
    values are 1D numpy arrays.
    """
    fig = plt.figure()
    ax = fig.gca()
    #axins = inset_axes( ax, width='70%', height=1., loc=10 )
    zoom = 2.5
    axins = zoomed_inset_axes( ax, zoom, 1 ) # location == upper right 

    # create a color instance
    rcolors = numpy.random.random( (len(cells),3) )
    cc = colors.ColorConverter()
    cmap = [ cc.to_rgb( c ) for c in rcolors ]

    nx_max = ny_max = 0

    for i, data in enumerate( ts.itervalues() ):
        thebins = numpy.arange( 0, data.max(), 4 )
        hist, bins = numpy.histogram( data, bins=thebins )
        # store for later use
        nx_max = max( thebins.max(), nx_max )
        if log:
            hist = numpy.log1p( hist )
        ny_max = max( hist.max(), ny_max )
        ax.plot( bins[:-1], hist, 'o-', ms=5, color=cmap[i] )
        axins.plot( bins[:-1], hist, 'o-', ms=6, color=cmap[i] )

    ax.set_xlabel( r'Generator lifespan', fontsize=20 )
    if log:
        ax.set_ylabel( r'Number of generators ($\log_{10}$)', fontsize=20 )
    else:
        ax.set_ylabel( r'Number of generators', fontsize=20 )

    # find axis max over all time series
    # ny_max = max( [ y.max() for y in ts.values() ] )
    # sub region of the original image
    yscale = ny_max * (0.4)
    # set the zoom cutoff as a percentage of the max
    xmax = cutoff * nx_max
    x1, x2, y1, y2 = 10, xmax, 0, yscale

    print x1, x2, y1, y2
    
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_aspect( xmax/yscale )
    axins.set_xticks([])
    axins.set_yticks([])
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    fig.show()
    return fig
 
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


def plot_midrange_ts( new_file, old_file, skip=10, fontsize=20,
                      lines=None, plot_mean=False ):
    """
    Plots two times series of midrange generators.
    """
    ts_new = numpy.loadtxt( new_file )
    ts_old = numpy.loadtxt( old_file )

    fig = plt.figure( dpi=160, figsize=([10,4]) )
    ax = fig.gca()
    ax.plot( ts_new[::skip], 'b-', linewidth=2 )
    ax.plot( ts_old[::skip], 'r-', linewidth=2 )
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

    # plot mean and std range
    if plot_mean:
        
        # for plotting average -- these are already log values if
        # log==True in kwargs
        yhist = numpy.array( all_ny, dtype=numpy.float64 )#.ravel()
        avg = yhist.mean( axis=0 )
        err = yhist.std( axis=0 )
        upper = avg + err
        lower = avg - err

       # now plot the interpolated average distribution here so it is on
        # top of the other stuff
        yp = None
        xi = numpy.linspace( 0, bins[-1],200)
        yi = stineman_interp( xi, bins[:-1], avg, yp )
        # interpolate upper and lower error bars to get envelope
        # COMMENTED OUT BELOW
        upper_yi = stineman_interp( xi, bins[:-1], upper, yp )
        lower_yi = stineman_interp( xi, bins[:-1], lower, yp )
    
    
    fig.show()
    return fig

def plot_hist_cut_axis( cells, color='blue',
                        normed=False, fontsize=20,
                        cell_type='New', skip=1,
                        show_plot=False, log=False,
                        left_xlim=300, right_xlim=1600,
                        cutoff=0.08, ts_max=100, histtype='bar' ):
    """
    Plot a histogram of generator lifespans computed as distance from diagonal.

    cell -- full path to perseus output file.
    """
    from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes
    from mpl_toolkits.axes_grid.inset_locator import mark_inset

    # timing
    start = time.time()
   
    # cell is a list of frames
    # create an array with first entry in cell, then extend in
    # loop
    # ts = [ get_ts( cell[0] ) ]
    # #    ts = ts[:-1]
    # for f in cell[1:ts_max]:
    #     ts.append( get_ts( f ) )
    #     #new_ts = get_ts ( f )
    #     # the (almost) infinite generator overwhelms the plot
    #     #      new_ts = new_ts[:-1]
    #     # extend last ts of generators by new_ts
    #     #ts = numpy.hstack( ( ts, new_ts ) )
    #     # the histogram of the data
    #     #vals = vals[1:]
    # # convert and normalize
    # ts = numpy.hstack( ts )
    # ts = numpy.asarray( ts, dtype=numpy.float )
    # ts /= ts.max()

    all_ts = []
    for c, cell in enumerate( cells ):
        # list of lifespans for each frame in cell
        # concatenate all of the diagrams for cell lifespans into one
        # timeseries
        ts = [ get_ts ( frame ) for frame in cell[:ts_max:skip] ]
        ts = numpy.hstack( ts )
        # convert and normalize, then append to the list of time series
        ts = numpy.asarray( ts, dtype=numpy.float )
        if not log:
            ts /= ts.max()
        all_ts.append( ts )

        newtime = time.time()
        print "Done creating time series. Took", newtime - start, " seconds"
        start = newtime

    # This method works with mpl 0.99 ( plt.subplot() works with >1.1.0 )
    fig = plt.figure( )
    ax = fig.add_subplot( 121 ) 
    ax2 = fig.add_subplot( 122 )

    # create two identical histograms 
    n, bins, patches = ax.hist( all_ts, bins=int(len(ts)/10.),\
                                normed=normed, facecolor=color,\
                                alpha=0.75, histtype=histtype, log=log )
    n, bins, patches = ax2.hist( all_ts, bins=int(len(ts)/10.),\
                                 normed=normed, facecolor=color,\
                                 alpha=0.75, histtype=histtype, log=log )

    # make the inset
    # options of loc: BEST, UR, UL, LL, LR, R, CL, CR, LC, UC, C = range(11)
    #axins = zoomed_inset_axes(ax, 10, loc=10 ) # ( axes, zoom power, location )
    if not log:
        axins = inset_axes( ax, width='70%', height=1., loc=10 )
        axins.hist( ts, bins=int(len(ts)/10.),\
                    normed=normed, facecolor=color,\
                    alpha=0.75, histtype=histtype, log=log )

    print "Done creating histograms. Took", time.time() - start, " seconds"
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
    ax.plot( (1-d,1+d),(-d,+d), **kwargs )      # bottom-right diagonal
    ax.plot( (1-d,1+d),(1-d,1+d), **kwargs )    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the right axes
    ax2.plot((-d,+d),(-d,+d), **kwargs)   # bottom-left diagonal
    ax2.plot((-d,+d),(1-d,1+d), **kwargs) # top-left diagonal

    # vertical separator
    ax2.set_ylabel( '---------- ---------- ---------- ---------- ----------', horizontalalignment='center' )

    x_max = max( [x.max() for x in all_ts] )
    y_max = max( [y.max() for y in n ] )
    # zoom-in / limit the view to different portions of the data
    ax.set_xlim( 0., left_xlim ) # most of the data
    #ax.set_ylim( 0., y_max+ 10 )
    ax2.set_xlim( right_xlim, x_max + 0.01) # outliers/inf gens
    #ax2.set_ylim( 0., 0.001 * y_max )  # zoom in to see the inf gens

    if not log:
        # sub region of the original image
        yscale = n.max() * (1./60)
        # set the zom cutoff as a percentage of the max
        xmax = cutoff * ts.max()
        x1, x2, y1, y2 = 0.015, xmax, 0, yscale
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_aspect( xmax/yscale )
        axins.set_xticks([])
        axins.set_yticks([])

    # set some tick marks
    left_ticks = numpy.arange( 0, left_xlim-0.01, 0.02 )
    left_labels = [ str( x ) for x in left_ticks ]
    ax.set_xticks( left_ticks )
    ax.set_xticklabels( left_labels )
    # scientific notation
    if not log:
        ax.ticklabel_format( style='sci', scilimits=(0,0), axis='y' )
    
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    if not log:
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    fig.subplots_adjust( wspace=0.1 )
    plt.show()
    
    return fig, ts

def plot_hist_stack(  cells, color='blue',
                      normed=False, fontsize=20,
                      cell_type='New',
                      show_plot=False, log=False,
                      left_xlim=300, right_xlim=1600,
                      cutoff=0.08, ts_max=100, histtype='stepfilled',
                      rwidth=1, skip=1, nbins=500, nbins_inset=200):
    """
    cells -- list of cells (full paths to) whose generator lifespans
    we want to stack in a single histogram.

    This is similar to the above, plot_hist_cut_axis(), except that it
    stacks numerous histograms.
    """
    from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes
    from mpl_toolkits.axes_grid.inset_locator import mark_inset
    import matplotlib.colors as colors


    # create a color instance
    rcolors = numpy.random.random( (len(cells),3) )
    cc = colors.ColorConverter()
    cmap = [ cc.to_rgb( c ) for c in rcolors ]

    # This method works with mpl 0.99 (plt.subplot() works with >1.1.0 )
    fig = plt.figure( figsize=(6,5) )
    ax = fig.add_subplot( 121 ) 
    ax2 = fig.add_subplot( 122 )

    # time shit
    t0 = start = time.time()

    # cell == cell directory containing persistence diagrams
    all_ts = []
    for c, cell in enumerate( cells ):
        # list of lifespans for each frame in cell
        # concatenate all of the diagrams for cell lifespans into one
        # timeseries
        ts = [ get_ts ( frame ) for frame in cell[:ts_max:skip] ]
        ts = numpy.hstack( ts )
        # convert and normalize, then append to the list of time series
        ts = numpy.asarray( ts, dtype=numpy.float )
        ts /= ts.max()
        all_ts.append( ts )

        newtime = time.time()
        print "Done creating time series. Took", newtime - start, " seconds"
        start = newtime
        
    min_len = min( [ len( t ) for t in all_ts ] )
    n, bins, patches = ax.hist( all_ts, bins=nbins,#min_len/float(2*ts_max),
                                histtype='barstacked', rwidth=rwidth, color=cmap,
                                log=log )
    n2, bins2, patches2 = ax2.hist( all_ts, bins=nbins,#min_len/float(2*ts_max),
                                    histtype='barstacked', rwidth=rwidth, color=cmap,
                                    log=log )
    if not log:
        axins = inset_axes( ax, width='70%', height=1., loc=10 )
        n3, bins2, patches3 = axins.hist( all_ts, bins=nbins_inset, 
                                          normed=normed,\
                                          alpha=0.75, histtype='barstacked', log=log, \
                                          rwidth=rwidth, color=cmap )


    # color the patches
    if not log:
        all_patches = [ patches, patches2, patches3 ]
    else:
        all_patches = [ patches, patches2 ] 
    for patch_set in all_patches:
        # each patch_set contain a set of stacked recatcngles. 
        for ci, thispatch in enumerate( patch_set ):
            for p in thispatch:
                p.set_color( cmap[ci] )
   

    # for P in [patches, patches2, patches3]:
    #     for patch in P:
            
         
    print "Done creating histograms. Total time: ", time.time() - t0, "seconds"
       
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
    ax.plot( (1-d,1+d),(-d,+d), **kwargs )      # bottom-right diagonal
    ax.plot( (1-d,1+d),(1-d,1+d), **kwargs )    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the right axes
    ax2.plot((-d,+d),(-d,+d), **kwargs)   # bottom-left diagonal
    ax2.plot((-d,+d),(1-d,1+d), **kwargs) # top-left diagonal

    # vertical separator
    ax2.set_ylabel( '---------- ---------- ---------- ---------- ----------', horizontalalignment='center' )

    # zoom-in / limit the view to different portions of the data
    # We key some of the limits off the max over all time series
    nx_max = max( [ x.max() for x in all_ts ] )
    ny_max = max( [ y.max() for y in n ] )
    left_ymax = sum( [ y.max() for y in n ] ) + 10

    print "nymax", ny_max
    
    ax.set_xlim( 0., left_xlim ) # most of the data
    ax.set_ylim( 0., left_ymax )
    ax2.set_xlim( right_xlim, nx_max + 0.01) # outliers/inf gens
    ax2.set_ylim( 0., 0.002 * ny_max )  # zoom in to see the inf gens

    if not log:
        # sub region of the original image
        yscale = ny_max * (1./20)
        # set the zoom cutoff as a percentage of the max
        xmax = cutoff * nx_max
        x1, x2, y1, y2 = 0.015, xmax, 0, yscale
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_aspect( xmax/yscale )
        axins.set_xticks([])
        axins.set_yticks([])

    # set some tick marks
    left_ticks = numpy.arange( 0, left_xlim-0.01, 0.02 )
    left_labels = [ str( x ) for x in left_ticks ]    
    ax.set_xticks( left_ticks )
    ax.set_xticklabels( left_labels )
    if not log:
        ax.ticklabel_format( style='sci', scilimits=(0,0), axis='y' )
    
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    if not log:
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    fig.subplots_adjust( wspace=0.1 )

    plt.show()

    return fig, all_ts

def load_birth_times( old_hist = 'data/old_birth_times.pkl',
                      new_hist = 'data/new_birth_times.pkl' ):
    with open( old_hist ) as fh:
        old_bt = pkl.load( fh )
    with open( new_hist ) as fh:
        new_bt = pkl.load( fh )
    return old_bt, new_bt


def load_all_bt( prefix='' ):
    BT = defaultdict( dict )
    for eps in [30,40,50]:
        old_ = prefix+'./bt_old_normed_eps'+str(eps)+'.pkl'
        new_ = prefix+'./bt_new_normed_eps'+str(eps)+'.pkl'
        old, new = load_birth_times( old_, new_ )
        BT[eps]['old'] = old
        BT[eps]['new'] = new
    return BT
         

if __name__ == "__main__":

    if 1:
        prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
        #prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
        newlist = ['new_10' , 'new_110125', 'new_130125', 'new_140125', 'new_3',
                   'new_4', 'new_40125', 'new_50125', 'new_6', 'new_60125', 'new_9']
        oldlist = ['old_100125', 'old_120125', 'old_15', 'old_2', 'old_4000', 'old_4001',
                   'old_5',  'old_50125',  'old_6',  'old_7',  'old_8',  'old_9',  'old_90125']

        cells = [ prefix + c + '/' for c in newlist ]
        #cells = [ prefix + c + '/' for c in oldlist ]
        
        frames = [ dir_list( c ) for c in cells ]

        dc = dict_of_cell_dirs( cells )
        #ts = concat_timeseries( dc, ts_max=-1, skip=10 )

        # find thresholds where first generators above lifespan <eps> are born
        bt = []
        for eps in [30, 40]:
            print "eps = ", eps
            print ""
            for cell in cells:
                bt.append( get_birth_times( cell, eps ) )
            # save to disk!
            with open( './data/bt_new_normed_eps'+str(eps)+'.pkl', 'w' ) as fh:
                pkl.dump( bt, fh )
    
    #fig = plot_scatter( ts, log=True )
    
    # fig, ts = plot_hist_stack( frames, left_xlim=0.2, right_xlim=0.6, normed=False,
    #                            cutoff=0.2, ts_max=1000, skip=20, log=True )

    # plot_hist_cut_axis( frames, normed=False, log=True,
    #                     left_xlim=0.2, right_xlim=0.6,
    #                     cutoff=0.2, ts_max=1000, skip=20 )

    if 0:
        # birth time dict 
        BT = defaultdict( dict )
        for eps in [30,40,50]:
            old_ = './data/old_bt_eps'+str(eps)+'.pkl'
            new_ = './data/new_bt_eps'+str(eps)+'.pkl'
            old, new = load_birth_times( old_, new_ )
            BT[30]['old'] = old
            BT[30]['new'] = new
    
    if 0:
        # find the old maxes. do these in sequence to avoid killing
        # the disk with thousands of minute search (sigh).
        old_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
        old_dirs= os.listdir( old_prefix )
        old_cells = [ old_prefix + c + slash for c in old_dirs ]
        print "finding maxima for old cells..."
        with timer.Timer():
            old_maxes = all_maxes( old_cells )
            with open( './data/old_maxes.pkl', 'w' ) as fh:
                pkl.dump( old_maxes, fh )

        # now find the new maxes
        new_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
        new_dirs= os.listdir( new_prefix )
        new_cells = [ new_prefix + c + slash for c in new_dirs ]
        print "finding maxima for new cells..."
        with timer.Timer():
            new_maxes = all_maxes( new_cells )
            with open( './data/new_maxes.pkl', 'w' ) as fh:
                pkl.dump( new_maxes, fh )

      

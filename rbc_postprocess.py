"""
    
Module for assorted post processing of RBCs (Perseus Output)
created by kel 5/23/2012

    
    New Cells: 110125, 130125, 140125, 40125, 50125
    Old Cells: 100125, 120125, 50125, 90125
"""

import numpy
import matplotlib
import matplotlib.pyplot as plt
import re
import os
import rbc_current as rc

def plot_diagram (persFile, lb=0, ub=2, out_type='bin',rmv='Y', dpi=80, fontsize=20):
    """
    Plot persistence diagram for data in persFile. If rmv=='Y', remove
    the infinite generator (though this doesn't seem to be working).
    """
    with open(persFile, 'r') as fh:
        s = fh.read()
    fh.close()
    s = s.split('\n')#seperate gens
    s.remove('') #remove blank lines
#    goodGens = rc . get_gens_sigma(persFile,lb,ub)
#    goodGens = rc . get_gens_bin (persFile)
    goodGens = rc . get_outlier_gens (persFile, lb, ub, 
                                      out_type, rmv, '')
    ggList = zip(*goodGens)
    x = []
    y = []
    maxLevel = str(int(s[-1].split(' ')[-1])+1)
    for i in xrange(len(s)):
        birth,death = s[i].split(' ')
        if int(birth) == -1:
            birth = maxLevel
        if int(death) == -1:
            death = maxLevel
        if not (int(birth),int(death)) in goodGens:
            x.append( int(birth) )
            y.append( int(death) )
    fig = plt.figure( dpi=dpi )
    ax = fig.gca()
    ax.scatter( x, y,c='b',marker='o',lw=.1)
    ax.scatter( ggList[0], ggList[-1],c='r',marker='o',lw=.1)
    line = [0, int(maxLevel)]
    ax.plot(line, line, 'g-')
    ax.set_xlim( [0, max( x )+10] )
    ax.set_ylim( [0, max( y )+20] )
    ax.set_xlabel( r'birth', fontsize=fontsize )
    ax.set_ylabel( r'death', fontsize=fontsize )
    xticks = [ int( tk ) for tk in ax.get_xticks() ]
    yticks = [ int( tk ) for tk in ax.get_yticks() ]
    ax.set_xticklabels( xticks, fontsize=fontsize )
    ax.set_yticklabels( yticks, fontsize=fontsize )
    fig.show()
    return fig

def plot_diagram_std (persFile, fontsize=16, scale=1,
                      color='b', show_fig=True, fig=None,
                      shape='o', plot_inf=False ):
    """
    persFile -- path to <perseus output>_*.txt, where * is the dimension.

    scale -- Factor to scale the birth/death times. 
    """
    # cast values as floats for division
    s = numpy.loadtxt( persFile, dtype=numpy.float, delimiter=' ' )
    births = s[:,0]
    deaths = s[:,1]

    # max death time
    maxd = deaths.max()
    
    # non-infinite gens
    normal_idx = numpy.where( deaths != -1 )[0]
    # infinite gens
    inf_idx = numpy.where( deaths == -1 )[0]
    
    # scale to match another persistence diagram (used for undoing a
    # scaling applied to Gaussians in RBC paper).
    if scale != 1:
        s /= scale
        # reset stuff
        births = s[:,0]
        deaths = s[:,1]
        # max death time
        maxd = deaths.max()
        inf_vec = (maxd + 1) * numpy.ones( len( inf_idx ) )
    else:
        inf_vec = (maxd + 1) * numpy.ones( len( inf_idx ) )

    print "Max death time ",  maxd
    
    if not fig:
        fig = plt.figure( ) #dpi=160 )
        fig.patch.set_alpha( 0.0 )
    ax = fig.gca()

    # plot the normal generators
    ax.plot( births[normal_idx], deaths[normal_idx], color+shape )

    # create diagonal
    diag = [0, maxd+1]
    ax.plot(diag, diag, 'g-')

    # plot 'em
    if plot_inf:
        ax.plot( births[inf_idx], inf_vec, 'ro' )

    # fix the left x-axis boundary at 0
    xticks = [ int( tk ) for tk in ax.get_xticks() ]
    yticks = [ int( tk ) for tk in ax.get_yticks() ]
    ax.set_xticklabels( xticks, fontsize=fontsize )
    ax.set_yticklabels( yticks, fontsize=fontsize )
    ax.set_xlim( left=0 )
    if show_fig:
        fig.show()

    # total number of persistence intervals
    print "Total number of persistence intervals", len( births ) 
    return fig

def plot_diagram_regions( persFile, lines=None, fontsize=16, zoom=False, scale=None, gauss=False ):
    """
    persFile -- path to perseus persistence diagram text file

    lines  -- list of ints indicating region-separating h/vlines.

    ** NOTE ** The plot attributes are specifically for cell new11.
    """
    # with open(persFile, 'r') as fh:
    #     s = fh.read()
    # fh.close()
    # s = s.split('\n')#seperate gens
    # s.remove('') #remove blank lines
    # x = []
    # y = []
    s = numpy.loadtxt( persFile, dtype=numpy.float, delimiter=' ' )
    maxLevel = s.max()
    # locate the infinite generators
    w = numpy.where( s == -1 )
    if scale:
        s = numpy.asarray( s, dtype=numpy.float )
        s /= scale
    s[ w ] = s.max()

    nx = s[:,0]
    ny = s[:,1]
    
    # for i in xrange(len(s)):
    #     birth,death = s[i].split(' ')
    #     if int(birth) == -1:
    #         birth = maxLevel
    #     if int(death) == -1:
    #         death = maxLevel
    #     x.append( int(birth) )
    #     y.append( int(death) )
    # now make the figure
    fig = plt.figure()# dpi=160, frameon=False )
    ax = fig.gca()
    ax.scatter( nx, ny,c='b',marker='o',lw=.1, s=50)
    diag = [0, maxLevel]
    ax.plot( diag, diag, 'g-')
    if not gauss:
        xticks = [0,500,1800] +lines
        yticks = [0,500,1800] +lines
        # xticks = [0,500,1000,1500,2000,2500] 
        # yticks = [0,500,1000,1500,2000,2500]
    else:
        ## THESE WEIRD VALUES ARE FOR A GAUSSIAN WITH NOISE AND SUBPEAK
        xticks = [0,5,10,15,20] + lines
        #xticks.pop( xticks.index(19) )
        yticks = [0,5,10,15, 20] + lines
        # xticks = [0,500,1500,2000,2500] + lines
        # yticks = [0,500,1500,2000,2500] + lines
    xticks.sort()
    yticks.sort()
    xticks_str = [ str( t ) for t in xticks ]
    yticks_str = [ str( t ) for t in yticks ]
    ax.set_xticks( xticks )
    ax.set_yticks( yticks )
    ax.set_xticklabels( xticks_str, fontsize=fontsize )
    ax.set_yticklabels( yticks_str, fontsize=fontsize )
    if lines:
        for line in lines:
            ax.hlines( line, 0, line, linestyles='dashed' ) 
            ax.vlines( line, line, s.max()+1, linestyles='dashed' )
    if zoom:
        ax.set_xlim( (lines[0]-200, lines[1]+200) )
        ax.set_ylim( (lines[0]-200, lines[1]+200) )
        ax.set_autoscale_on( False )
    else:
        print max(nx)
        ax.set_xlim( [0, max( nx )+1] )
        ax.set_ylim( [0, max( ny )+1] )
    #ax.set_title( 'Persistence Diagram', fontsize=fontsize+4 )
    #ax.set_xlabel( 'birth', fontsize=fontsize )
    #ax.set_ylabel( 'death', fontsize=fontsize )
    fig.show()
    return fig

    

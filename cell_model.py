"""
Module to create a model of a RBC by approxiamting it using a low-band pass. WE take only the first N modes (N <= 3 probably).

Opened: Aug. 2, 2012

Author: Jesse Berwald
"""
import numpy as np
import pylab
from scipy import fftpack


def load_image( fname ):
    """
    """
    return np.loadtxt( fname )

def fft_image( img ):
    """
    Perform 2D fft on a cell frame and return the shifted (centered)
    version.

    img -- 2D array of floats ( eg. f(x,y) = z )
    """
    transform = fftpack.fft2( img )
    # shifts spatial frequencies to center of image
    return fftpack.fftshift( transform )

def ifft_image( X ):
    """
    Return to the spatial domain.
    """
    Y = fftpack.ifft2( X )
    return fftpack.ifftshift( Y )

def spectral_power( X, take_log=False ):
    if take_log:
        return np.log1p( np.abs( X ) )
    else:
        return np.abs( X )

def band_pass( X, dx, dy ):
    """
    Create band pass matrix H composed of the the frequencies in X in
    the squaree (center+dx), (center+dy)
    """
    nx,ny = X.shape
    # the largest value the index can attain
    max_idx = min( nx, ny )
    # locate the (row,col) center NOTE: l0 == "el zero"
    k0,l0 = nx/2, ny/2
    H = np.zeros_like( X )

    # Only keep the DC component
    if dx == 0 and dy == 0:
        H[k0,l0] = X[k0,l0]
        return H
    # else
    for i in xrange( k0-dx, k0+dx ):
        # now we walk down a column j from the ith row
        for j in xrange( l0-dy, l0+dy ):
            H[i,j] = X[i,j]
    return H


if __name__ == "__main__":

    from local_settings import *
  
    img = load_image( frame )
    orig = np.load( orig_frame )

    # slice row
    row = 100
    #num_modes = 4

    for i,num_modes in enumerate([4]):
        # fft for non-halo cell
        X = fft_image( img )
        X = spectral_power( X )
        H = band_pass( X, num_modes, num_modes )
        Y = ifft_image( H )
        Y = spectral_power( Y )

        bnd = np.loadtxt( bfile )

        # fft for original cell with halo
        Xhalo = fft_image( orig )
        Xhalo = spectral_power( Xhalo )
        Hhalo = band_pass( Xhalo, num_modes, num_modes )
        Yhalo = ifft_image( Hhalo )
        Yhalo = spectral_power( Yhalo )

        if 0:
            # draw stuff 
            fig1 = pylab.figure(1, dpi=160) # i, figsize=(10,8) )
            # boundary removed
            ax1 = fig1.add_subplot( 111 )
            ax1.imshow( img )
            ax1.set_xticks( [] )
            ax1.set_yticks( [] )
            # label the row
            # ax1.set_ylabel( 'Boundary removed', fontsize=18 )
            #ax1.set_title( 'Intensity function', fontsize=18 )
            # fig.colorbar( ax=ax1, cax=ax1 )
            fig1.show()

        if 0:
            # first few modes
            fig2 = pylab.figure(2, dpi=160)
            ax2 = fig2.add_subplot( 111 )
            ax2.imshow( Y )
            #ax2.set_title( r''+str( num_modes ) + ' modes' )
            ax2.set_xticks( [] )
            ax2.set_yticks( [] )
            fig2.show()

        if 1:
            # slices of img and Y
            height = 1400
            fig3 = pylab.figure(3, dpi=160)
            ax3 = fig3.add_subplot( 111 )
            ax3.plot( img[row,:], 'b--', lw=3, label=r'Intensity $f(x,y)$' )
            ax3.plot( Y[row,:], 'g-', lw=4, label=r'Approx.' )
            ax3.hlines( height-1, 1, 197, linestyle='dotted', linewidth=2 )        

            # find sublevel set
            # mask regions above height and fill those below
            y1 = Y[row,:]
            w = np.where( y1 <= height )[0]
            y2 = y1[w]
            y_masked = np.ma.masked_greater( y1, height )
            ax3.fill_between( range(0,Y.shape[1]), 0, y_masked, where=y_masked<=height,
                              facecolor='green', alpha=0.5, interpolate=True,
                              label=r'Sublevel set')
                #ax3.set_title( r'Horizontal cross-section', fontsize=24 )
            #        ax3.set_aspect( 'equal', 'box' )
            ax3.set_xlim( [0,197] )
            xticks = [ str( int(x) ) for x in ax3.get_xticks() ]
            yticks = [ str( int(y) ) for y in ax3.get_yticks() ]
            ax3.set_xticklabels( xticks, fontsize=20 )
            ax3.set_yticklabels( yticks, fontsize=20 )        
            leg = ax3.legend()
            for t in leg.get_texts():
                t.set_fontsize('small')
            fig3.show()


        # ## HALO ##
        # # original without boudnary removed
        # ax4 = fig.add_subplot( 234 )
        # ax4.imshow( orig )
        # ax4.set_ylabel( 'Orignal with Halo', fontsize=18 )

        # # approx cell using <halo_modes>
        # ax5 = fig.add_subplot( 235 )
        # ax5.imshow( Yhalo )

        # # slices
        # ax6 = fig.add_subplot( 236 )
        # ax6.plot( orig[row,:], 'b-', lw=2, label='cell (halo)' )
        # ax6.plot( Yhalo[row,:], 'r-', lw=2, label='FFT approx.' )
        # #ax6.set_title( r'Sliced along row='+str(row) )

        #fig.savefig( theDir + 'cell_approx_' + str(num_modes) + 'modes.png' )

        # fig1.show()
    

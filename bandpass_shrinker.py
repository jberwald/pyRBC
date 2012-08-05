import cPickle as pkl
import numpy as np

def run( fname ):
    with open( fname ) as fh:
        data = pkl.load( fh )
    mats = [ (frame['ifft_nonzero'], frame['mean']) for frame in data ]
    for i in range( len(data) ):
        data[i]['ifft_nonzero'] = np.ma.masked_less( mats[i][0], mats[i][1] ).mask
    # save back to disk with most efficient protocol
    with open( fname, 'w' ) as fh:
        pkl.dump( data, fh, protocol=-1 )
    print "Done shrinking", fname
        
if __name__ == "__main__":

    low = '00'
    high = '01'
    modes = np.linspace( 0, 1, 11 )
    for i,mode in enumerate( modes[:-1] ):
        low=str(mode)
        low_ = ''.join(low.split('.'))
        high=str(modes[i+1])
        high_ = ''.join(high.split('.'))
        filename = 'fft_r'+low_+'_r'+high_+'.pkl'
        print "shrinking", filename
        run( filename )

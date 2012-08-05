from numpy import *
from pylab import matshow


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
        -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def dist( x1, x2, y1, y2 ):
    return (x1-x2)**2 + (y1-y2)**2

def noisy_gaussian(height, center_x, center_y, width_x, width_y, noise_level, shape):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
        -(((center_x-x)/width_x)**2 + ((center_y-y)/width_y)**2)/2) + noise_level*random.random( shape )
       
def plot_gaussian( size_x, size_y, height, width=20, noise_level=None ):
    """
    Gaussian centered at (100,100). Equal width in each direction.
    """
    Xin, Yin = mgrid[0:size_x, 0:size_y]
    if not noise_level:
        data = gaussian(height, 100, 100, width, width)(Xin, Yin)
    else:
        data = noisy_gaussian(height, 100, 100, width, width, noise_level, shape=Xin.shape)(Xin, Yin) 
    zd = zeros_like( data )
    for i in range( data.shape[0] ):
        for j in range( data.shape[1] ):
            if dist( 100, i, 100, j ) < 10*width**2:
                zd[i,j] = 1
    data = zd * data
    if noise_level:
        save( "./data/gauss_bump_noise"+str(noise_level)+"_"+str(height), data )
    else:
        save( "./data/gauss_bump_"+str(height), data )
    matshow( data )
        
def make_gaussian( rows, cols, height=2, width_x=5, width_y=5, noise_level=None ):
    """
    Make an index array of size (rows,cols), and pass that into an
    instance of gaussian() or noisy_gaussian(). Returns a numpy array
    of size (rows,cols).

    Example: make_gaussian( 101, 101 ) returns a smooth gaussian
    centered at (50,50) with default width of 5 in the x and y
    directions and a height of 2 at the center.
    """
    idx = indices( (rows,cols) )
    if noise_level:
        f = noisy_gaussian( height,
                            int( rows/2. ),
                            int( cols/2. ),
                            width_x,
                            width_y,
                            noise_level=noise_level )
    else:
        f = gaussian(  height,
                       int( rows/2. ),
                       int( cols/2. ),
                       width_x,
                       width_y )
    gauss = f( idx[0], idx[1] )
    return gauss
    

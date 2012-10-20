import numpy as np
import matplotlib.pyplot as plt

def perturb_function( nx, vf, nsamp=5 ):
    """
    Given f, perturb a few points, but keep the two functions close.

    nx -- x points

    vf -- vectorized version of f
    """
    # get indices to perturb
    sample_idx = np.random.random_integers( 0, len(nx), size=(nsamp,) )

    print sample_idx

    ny = vf( nx )
    # perturb ny
    ny[ sample_idx ] += np.random.sample( size=(nsamp,) )

    return ny

def plot_func( nx, ny, fig=None, diag=True ):
    """
    """
    if not fig:
        fig = plt.figure()
        color = "b"
        label = r"$f(x)$"
    else:
        print "hello"
        color = "r"
        label = r"$g(x)$"
    if not diag:
        ax = fig.gca()
        ax.plot( nx, ny, lw=3, color=color, label=label )
    else:
        ax1 = fig.add_subplot( 121 )
        ax1.plot( nx, ny, lw=3, color=color, label=label )
        ax1.set_xlabel( r'$x$', fontsize=20 )
        ax1.set_ylabel( r'$f(x)$', fontsize=20 )
        # plot diagonal
        ax2 = fig.add_subplot( 122 )
        ax2.plot( nx, nx, 'k-', lw=1.5 )
        ax2.set_xlabel( r'$\tau$', fontsize=20 )
    fig.show()
    return fig

def main():
    # create a couple of polynomials that are "close".
    nx = np.linspace( 0, 5, 200 )
    f1 = lambda x : 0.3 * ( x * (0.9*x - 1) * (0.8*x - 2.5) * (x - 4) * (x-4.5) ) + 3
    f2 = lambda x : 0.3 * ( x * (0.9*x - 1) * (0.8*x - 2.5) * (x - 4) * (x-4.5)*\
                            (0.9*x - 0.2) * (0.5*x -2) * (0.7*x - 2.1) ) + 3
    vf1 = np.vectorize( f1 )
    vf2 = np.vectorize( f2 )
    ny1 = vf1( nx )
    ny2 = vf2( nx )
    thefig = plot_func( nx, ny1 )
    #thefig = plot_func( nx, ny2, fig=thefig )

    plt.legend()

    return thefig 


if __name__ == "__main__":

    fig = main()

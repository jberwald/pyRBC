"""
Broken axis example, where the y-axis will have a portion cut out.
"""
import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import ScalarFormatter

# 30 points between 0 0.2] originally made using np.random.rand(30)*.2
# pts = np.array([ 0.015,  0.166,  0.133,  0.159,  0.041,  0.024,  0.195,
#     0.039, 0.161,  0.018,  0.143,  0.056,  0.125,  0.096,  0.094, 0.051,
#     0.043,  0.021,  0.138,  0.075,  0.109,  0.195,  0.05 , 0.074, 0.079,
#     0.155,  0.02 ,  0.01 ,  0.061,  0.008])

pts = np.random.rand(30) 

# Now let's make two outlier points which are far away from everything.
#pts[[3,14]] += .8

# If we were to simply plot pts, we'd lose most of the interesting
# details due to the outliers. So let's 'break' or 'cut-out' the y-axis
# into two portions - use the top (ax) for the outliers, and the bottom
# (ax2) for the details of the majority of our data
#f,(ax,ax2) = plt.subplots(2,1,sharex=True)

fig = plt.figure()
ax = fig.add_subplot( 121 )
ax2 = fig.add_subplot( 122, sharey=ax )

# plot the same data on both axes
# ax.plot(pts)
# ax2.plot(pts)
ax.hist( pts, bins=10 )
ax2.hist( pts, bins=10 )

# make the inset
# options of loc: BEST, UR, UL, LL, LR, R, CL, CR, LC, UC, C = range(11)
axins = inset_axes( ax, width='50%', height=1., loc = 10 )
#axins = zoomed_inset_axes(ax, 1.0, loc=10) # zoom = 6
axins.hist( pts )

# zoom-in / limit the view to different portions of the data
ax2.set_xlim(.78,1.) # outliers only
ax.set_xlim(0,.22) # most of the data
ax.set_ylim(0,10)

# sub region of the original image
x1, x2, y1, y2 = 0, 0.15, 0, 3
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# hide the spines between ax and ax2
ax.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax.xaxis.tick_bottom()
ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d,+d),(-d,+d), **kwargs)      # top-left diagonal
ax.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-right diagonal

ax.ticklabel_format( style='sci', axis='y' )  # scientific notation

sf = ScalarFormatter()
sf.set_scientific( True )

# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

fig.subplots_adjust( wspace=0.1 )

plt.show()


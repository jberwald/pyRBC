import numpy as n
from pylab import figure, show
import matplotlib.cm as cm
import matplotlib.colors as colors

fig = figure()
ax = fig.add_subplot(111)
Ntotal = 1000
N, bins, patches = ax.hist(n.random.random((Ntotal,2)), 20)

#I'll color code by height, but you could use any scalar

# we need to normalize the data to 0..1 for the full
# range of the colormap
fracs = [n.astype(float)/n.max() for n in N]
norm = [colors.normalize(frac.min(), frac.max()) for frac in fracs ]

for i, (thisfrac, thispatch) in enumerate( zip(fracs, patches) ):
    color = cm.jet(norm[i](thisfrac))
    for p in thispatch:
        #thispatch.set_facecolor(color)
        p.set_facecolor(color)

show()

#!/usr/bin/python
import sys, os, subprocess
from subprocess import call
import rbc_npy2Perseus as n2p
"""
Make a movie from PNG files. Assumed that this script resides in the
folder containing a sequence of frames (PNG files).
"""

slash = "/"

def sort_files( dlist ):
    """ Sort files in dlist according to the order
    of creation.  """
    return sorted(dlist, key=order)

def order( line ):
    """split line and return the index of the partial cover"""
    return int( line.rstrip().split('.')[-2] )

def prefix_files( dlist, prefix ):
    """Return only files beginning with <prefix>"""
    plist = [x for x in dlist if x.startswith( prefix )]
    return plist

def suffix_files( dlist, suffix ):
    """Return only files ending with <suffix>"""
    slist = [x for x in dlist if x.endswith( suffix )]
    return slist

def patch_cover( dlist, prefix, suffix='png'):
    """Return image (of type <suffix>) that begin with <prefix>"""
    a = prefix_files( dlist, prefix )
    return suffix_files( a, suffix ) 

def run(dir, output, prgm='mencoder'):
    if not output.endswith('.mp4'):
        output+='.mp4'
    if not dir.endswith('/'):
        dir+='/'
    if os.path.isdir (dir):
        dlist = os.listdir(dir)
        frames = []
        for f in dlist:
            if f.endswith('npy') and not os.path.isdir(dir+f):
                frames.append(dir+f)
    else:
        print 'Error - input is not a directory'
        return
    #print frames
    frames.sort(key=n2p.natural_key)
    print '..creating movie..'
    if prgm.startswith('m'):
        call(['mencoder',
           frames,
           '-mf',
           'type=png:h=800:w=400:fps=34',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=avi',
           '-oac',
           'copy',
           '-o',
           output])
    else:
        call(['ffmpeg', '-qscale 5', '-r 34', '-b 9600', 'data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/new_110125-concatenated-ASCII_%.npy', output ])
    print '..finished..'

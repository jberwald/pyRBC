#!/usr/bin/env python
"""

Run chomp on .npys
"""

import numpy
import Image
import rbc_basic as R
import fft_image as FFT
import chomp_betti as CB
import os
import matplotlib.pyplot as plt
#import convert2vtk as cc
from pylab import *

#fname = '/data/jberwald/wyss/data/Cells_Jesse/New/cells/new_140125-concatenated-ASCII'
#file = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_140125/new_140125-concatenated-ASCII_24.npy'
#path = '/data/jberwald/wyss/data/Cells_Jesse/frames/new_140125/fft_frames/'
#bnd_file = '/data/jberwald/wyss/data/Cells_Jesse/New/cells/boundary_Nov_new140125'
#file = '/home/kellys/Dropbox/Shared-MSU-Kel/Research/RBC/frames/new_140125-concatenated-ASCII_24.npy'
#t_file = '/home/kellys/Dropbox/Shared-MSU-Kel/Research/RBC/frames/fft_frames/normed_frames/r05_cub/new_140125-concatenated-ASCII_0_r05.cub'

#Extracting fft code
#files = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/fft_frames/r10/new_110125-concatenated-ASCII_99_r10.pkl'
#files = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/fft_frames/r01/'
#files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/fft_frames/r01/'
#files = '/home/kellys/RBC_Deriv_120125_r01/'
#files = '/home/kellys/RBC_Deriv_110125_r01/'
#files = '/home/kellys/RBC_Deriv_120125_r005/'
#files = '/data/rbc_deriv_files/RBC_Deriv_120125_r01/'
#subdir = 'betti/'

#files = '/sciclone/home04/jberwald/data10/wyss/data/Cells_Jesse/New/frames/new_40125/'\
#   'fft_frames/bandpass/smooth_r0708/'
files = '/sciclone/home04/jberwald/data10/wyss/data/Cells_Jesse/Old/frames/old_50125/'\
    'fft_frames/bandpass/smooth_r0708/'

#grab files, skip directories
if os.path.isdir (files):
	fdir = files + '/'
	dlist = os.listdir(fdir)
	frames = []
	for f in dlist:
		if f.endswith('npy') and not os.path.isdir(fdir+f):
			frames.append(f)
else:
	frames = [ files ]
	fdir = files.rpartition('/')[0].rpartition('/')[0]+'/'
#savedir = '/home/kellys/RBC_Deriv_120125_r01/'
#savedir = '/home/kellys/RBC_Deriv_110125_r01/'
savedir = files + 'betti/' #'/home/kellys/RBC_Deriv_120125_r005/'

print 'attempting chomp runs..'
for frame in frames:
	filename = frame.rstrip('npy') + 'cub'
	savename = frame.rstrip('npy') + 'cbetti'
	#CB.png2chomp(fdir + frame)
	arr = numpy.load(fdir+frame)
	w = numpy.where(arr==True)
	#w2 = numpy.where(w[2] == 0)[0]
	newarr = numpy.vstack( (w[0],w[1])).T
	CB.array2chomp(newarr,fdir+filename)
	CB.run_chomp(fdir + filename, fdir+savename)
print 'finished..'


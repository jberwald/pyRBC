#!/usr/bin/env python
"""
Create derivative images

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
files = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_120125/fft_frames/r005/'
#files = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/fft_frames/r005/'

#grab files, skip directories
if os.path.isdir (files):
	fdir = files + '/'
	dlist = os.listdir(fdir)
	frames = []
	for f in dlist:
		if f.endswith('pkl') and not os.path.isdir(fdir+f):
			frames.append(f)
else:
	frames = [ files ]
	fdir = files.rpartition('/')[0].rpartition('/')[0]+'/'
#savedir = '/home/kellys/RBC_Deriv_120125_r01/'
#savedir = '/home/kellys/RBC_Deriv_110125_r01/'
savedir = '/home/kellys/RBC_Deriv_120125_r005/'

print 'attempting array creation..'
k= 0;
derivArr = [];
parArr = [];
frames.sort(key=R.natural_key)
for frame in frames:
	k = k+1
	savename = frame.rstrip('.pkl')
	D = FFT.extract_ifft(fdir+frame)
	arr = D['ifft_nonzero']
	masked_ifft = numpy.ma.masked_less_equal(arr,D['mean'])
	if(k != 1):
		savename = frame.rstrip('.pkl')
		newFile = numpy.abs(masked_ifft - prevFrame)
		m = newFile.mean()
		newFile = numpy.ma.masked_less_equal(newFile, m).mask
		numpy.save(savedir+savename + '.npy', newFile)
		#plt.imshow(newFile)
		#plt.savefig(savedir+savename+'.png')
		print k
		#derivArr.append(math.fabs(prevFrame[50][100] - masked_ifft[50][100]))
		#parArr . append(k)
	prevFrame = masked_ifft
	prevMean = D['mean'] #not currently used
	#print k


import re
def natural_key(string_):
    #Use with frames.sort(key=natural_key)
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]



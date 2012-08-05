#!/usr/bin/env python
"""
Module for creating movies from thresholded FFT, as well as analyzing first difference derivatives

Notations:

files is typically the path to the cell directory on simplex

savedir is new directory to save to

create_diff_image creates the first difference images for each frame in a cell
Note, on the labeling scheme: for the first

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

def create_diff_image ( files , savedir ):
	""" grab files, skip directories, unabashedly stolen from Jesse's code
        e.g. new cell name on simplex:
        '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/new_110125-concatenated-ASCII_324.npy'
        
    """
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
	
	#print 'attempting array creation..'
	k= 0;
	n=5000; #grab a selection of frames
	derivArr = [];
	parArr = [];
	frames.sort(key=natural_key)
	for frame in frames:
		k = k+1
		savename = frame.rstrip('.pkl')
		D = FFT.extract_ifft(fdir+frame)
		arr = D['ifft_nonzero']
		masked_ifft = numpy.ma.masked_less_equal(arr,D['mean'])
		if(k != 1 and k < n):
			savename = frame.rstrip('.pkl')
			newFile = numpy.abs(masked_ifft - prevFrame)
			m = newFile.mean()
			newFile = numpy.ma.masked_less_equal(newFile, m).mask
			numpy.save(savedir+savename + '.npy', newFile)
			#print k
		prevFrame = masked_ifft


def diff_pkl2chomp ( files , savedir ) :
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

	#print 'attempting chomp runs..'
	frames.sort(key=natural_key)
	for frame in frames:
		savename = frame.rstrip('npy') + 'cub'
		#CB.png2chomp(fdir + frame)
		arr = numpy.load(fdir+frame)
		w = numpy.where(arr==True)
		newarr = numpy.vstack( (w[0],w[1])).T
		CB.array2chomp(newarr,fdir+savename)
		CB.run_chomp(fdir + savename, savedir+savename)
	#print 'finished..'

def save_masks ( files , savedir ):
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
	#savedir = '/home/kellys/RBC_120125_old_r005/'

	frames.sort(key=natural_key)
	for frame in frames:
		savename = frame.rstrip('.pkl')
		#print 'preparing to extract..'
		D = FFT.extract_ifft( fdir+frame )
		arr = D['ifft_nonzero']
		masked_ifft = numpy.ma.masked_less_equal(arr, D['mean'] )
		plt.imshow(masked_ifft.mask)
		plt.savefig(savedir+savename+'.png')

import re
def natural_key(string_):
    #Use with frames.sort(key=natural_key)
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
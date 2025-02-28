###################################################################################
##########	For plotting CMDs, along with mass tracks 		########### 
###################################################################################

import re
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as img
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from scipy.interpolate import interp1d
from astropy.io.votable import parse
from numpy import nanmedian, nanstd, nanmean, median, std, mean, log10
import pandas as pd
import os

cd = os.chdir
pwd = os.getcwd
pd.options.mode.chained_assignment = None

from XRBID.DataFrameMod import Find, BuildFrame
from XRBID.Sources import LoadSources
from XRBID.Headers import heads, B, V, I, U, BV, VI, BI, UB, Filter, ID, X, Y
#from XRBID.DataFrameMod import SaveDF

default_aps = [0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.]

#min1 = pd.read_csv("masstrackmin_1Msun.frame", verbose=False)
#min2 = pd.read_csv("masstrackmin_2Msun.frame", verbose=False)
#min3 = pd.read_csv("masstrackmin_3Msun.frame", verbose=False)
#min5 = pd.read_csv("masstrackmin_5Msun.frame", verbose=False)
#min8 = pd.read_csv("masstrackmin_8Msun.frame", verbose=False)
#min20 = pd.read_csv("masstrackmin_20Msun.frame", verbose=False)
#min40 = pd.read_csv("masstrackmin_40Msun.frame", verbose=False)

#min_masses = [min1, min2, min3, min5, min8, min20, min40]
#min_masses = [min3, min5, min8, min20]

# Labels associated with each line. Sources BELOW each line (but above the previous) get this label
mass_labels = ["Low", "Intermediate", "Intermediate", "High"]

#print("Testing load on CMDs.py")


###-----------------------------------------------------------------------------------------------------

def MakeCMD(sources=None, xcolor=None, ycolor=None, xmodel=None, ymodel=None, figsize=(6,4), xlim=None, ylim=None, color="black", size=10, marker=None, label=None, save=False, savefile=None, title=None, subimg=None, annotation=None, annotation_size=None, imshow=True, fontsize=15, shift_labels=[[0,0],[0,0],[0,0],[0,0],[0,0]], set_labels=None, instrument="ACS", color_correction=[0,0], labelpoints=False, file_dir="/home/qiana/Documents/Research/"): 

	"""Makes a CMD from a given set of points, either from a list or an input dataframe.

	PARAMETERS: 
	sources [pd.dataframe or list] (None): 	Input may either be a pandas dataframe containing the appropriate magnitudes required for the CMD, 
				       		or a list of coordinates in the format [[xs],[ys]]. 
	xcolor [str or list]: 			The name of the color or magnitude to plot in the x-axis, as it is labeled in the dataframe. This will be the default x-axis label. If a list is given, it is 							assumed to contain the column names of the input dataframe to be subtracted (e.g. F555W - F814W)
	ycolor [str or list]: 			The name of the color or magnitude to plot in the y-axis, as it is labeled in the dataframe. This will be the default y-axis label. If a list is given, it is 							assumed to contain the column names of the input dataframe to be subtracted (e.g. F555W - F814W)
	xmodel [str or list]: 			The magnitude(s) of the filter(s) to be used from the stellar models for the x-axis. If given as a list, it is assumed the color is xmodel[0] - xmodel[1] 
	ymodel [str or list]: 			The magnitude(s) of the filter(s) to be used from the stellar models for the y-axis. If given as a list, it is assumed the color is ymodel[0] - ymodel[1] 
	figsize [tuple] (6,4): 			The desired dimensions of the figure. 
	xlim [tuple] (None): 			The limits on the x-axis. If none are given, the limits are assumed to be the (xmin - 1, xmax + 1) 
	ylim [tuple] (None) 			The limits on the y-axis. If none are given, the limits are assumed to be (ymax + 1, ymin - 1)
	color [str] ("black"): 			Marker color
	size [int] (10): 			Marker size
	marker [str] (None):			The style of the marker. Defaults to a filled point. If "o" is given, the marker is set to be an open circle. 
	label [str] (None):			The legend label to assign to the input points from sources. 
	save (bool) (False): 			Sets whether to automatically same the CMD image. 
	savefile [str] (None): 			The name to assigned the saved CMD image. 
	title [str)] (None): 			Title of the figure, to be placed near, but not at, the top of the figure.
	subimg [str] (None): 			Filename of an image to include in the corner of the CMD. This is to allow a subplot of the XRB plotted to be shown within the CMD. 
	annotation [str] (None): 		Additional annotation to add to the bottom corner of the CMD (usually XRB ID)
	annotation_size [int] (None): 		Annotation fontsize 
	imshow [bool] (True): 				
	fontsize [int] (20): 			The fontsize of text other than the annoations
	shift_labels [list]: 			List of x and y coordinate distances by which to shift each of the model mass labels.
						Defaults to [[0,0],[0,0],[0,0],[0,0],[0,0]] (no shifts)
	set_labels [list] (None):		Sets the position of labels. If none is given, positions are automatically calculated.
	instrument [str] ("ACS"): 		Name of the instrument used, to determine which models to call. 
	color_correction [list] ([0,0]): 	Corrections on the x and y position of the sources. Defaults to no correction. 
	labelpoints [list] (False): 		Labels to add to each point. If none are given, defaults to False and no labels added. 
	file_dir [str]: 			The directory within which the models may be found.

	RETURNS: 
	f, ax: 		Arguments defining the figure, which can be used to add more points to the CMD after the initial plotting.
 
"""

	# Setting the style of my plots
	#fontparams = {'font.family':'stix'}
	#labelparams = {'family':'stix', 'size':fontsize}

	temp_dir = pwd()
	cd(file_dir)

	# Reading in the appropriate models based on the instrument given.
	if instrument.upper() =="WFC3":
		mass1 = pd.read_csv("isoWFC3_1Msun.frame")
		mass3 = pd.read_csv("isoWFC3_3Msun.frame")
		mass5 = pd.read_csv("isoWFC3_5Msun.frame")
		mass8 = pd.read_csv("isoWFC3_8Msun.frame")
		mass20 = pd.read_csv("isoWFC3_20Msun.frame")
	elif instrument.upper() =="ACS": 
		mass1 = pd.read_csv("isoACS_WFC_1Msun.frame")
		mass3 = pd.read_csv("isoACS_WFC_3Msun.frame")
		mass5 = pd.read_csv("isoACS_WFC_5Msun.frame")
		mass8 = pd.read_csv("isoACS_WFC_8Msun.frame")
		mass20 = pd.read_csv("isoACS_WFC_20Msun.frame")

	masses = [mass1, mass3, mass5, mass8, mass20] # list of DataFrames of each mass model
	mass_labels = ["1 M$_\odot$", "3 M$_\odot$", "5 M$_\odot$", "8 M$_\odot$", "20 M$_\odot$"]

	if savefile: save = True

	# Setting the x- and y-axes labels.
	# if xcolor and/or ycolor is a list, will need to set color[0] - color[1] as the color of the appropriate axis
	if not xcolor: xcolor=xmodel
	if not ycolor: ycolor=ymodel

	if isinstance(xcolor, list): xlabel = " - ".join(xcolor)
	else: xlabel = xcolor
	if isinstance(ycolor, list): ylabel = " - ".join(ycolor)
	else: ylabel = ycolor

	### Pulling the x and y values of the sources ###
	# if input source is a pandas dataframe, read in the appropriate colors and magnitudes (and add correction, if needed)
	if isinstance(sources, pd.DataFrame):
		if isinstance(xcolor, list): xsources = sources[xcolor[0]].values - sources[xcolor[1]].values + color_correction[0]
		else: xsources = sources[xcolor].values + color_correction[0]
		if isinstance(ycolor, list): ysources = sources[ycolor[0]].values - sources[ycolor[1]].values + color_correction[1]
		else: ysources = sources[ycolor].values + color_correction[1]
		
	else: # If sources is a list or coordinates, pull the x and y values as given (with additional color correction)
		xsources = (np.array(sources[0]) + color_correction[0]).tolist()
		ysources = (np.array(sources[1]) + color_correction[1]).tolist()	
	### Will only need to call xsources or ysources from now on ###

	
	### PLOTTING MODEL MASS TRACKS ###

	f, ax = plt.subplots(figsize=figsize)

	# Setting the tick parameters
	ax.tick_params(direction="in", labelsize=15, bottom=True, \
		       top=True, left=True, right=True, length=7, width=2)
	ax.grid(alpha=0.8, linestyle="--")

	# If no xmodel or ymodel are given, default to that of the sources. 
	# This point will fail if the filters of the sources are not given a name corresponding to the filters in the model DataFrame!
	if not xmodel: xmodel=xcolor
	if not ymodel: ymodel=ycolor

	# Pulling the correct colors/magnitudes from the models and plotting
	# and setting the default mass track label position and plot limits
	xlims = []
	ylims = []
	for m, mass in enumerate(masses): # for each of the mass models, pull the color/mag given by xmodel and ymodel
		if isinstance(xmodel, list): 
			xtemp = mass[xmodel[0]].values - mass[xmodel[1]].values
			xtemp_label = max(xtemp) + 0.1 + shift_labels[m][0] # default mass label position, unless set_labels is given
			xtemp_left = max(xtemp)	# Keeping track of the leftmost x coordinate		
			if not xlim: xlims.append([min(xtemp)-1, max(xtemp)+1])
			invert_xlim = False
		else: # if x-axis is a magnitude..
			xtemp = mass[xmodel].values
			xtemp_label = min(xtemp) - 0.1 + shift_labels[m][0]
			xtemp_left = min(xtemp)
			if not xlim: xlims.append([max(xtemp)+1, min(xtemp)-1])
			invert_xlim = True
		if isinstance(ymodel, list): 
			ytemp = mass[ymodel[0]].values - mass[ymodel[1]].values
			if not ylim: ylims.append([min(ytemp)-1, max(ytemp)+1])
			invert_ylim = False
		else: # of y-axis is a magnitude...
			ytemp = mass[ymodel].values
			if not ylim: ylims.append([max(ytemp)+1, min(ytemp)-1])
			invert_ylim = True

		# Finding the best y-coordinate for the model label based on the leftmost
		ytemp_label = ytemp[xtemp.tolist().index(xtemp_left)] + invert_ylim*0.5 + shift_labels[m][1]

		# If set_labels is given, use this as the coordinate of the label
		# (overrides the label positions set above)
		if set_labels: 
			xtemp_label = set_labels[m][0]
			ytemp_label = set_labels[m][1]

		# Plotting mass track and mass label
		plt.plot(xtemp, ytemp, color="black", lw=1)
		ax.annotate(mass_labels[m], xy=(xtemp_label, ytemp_label), size=15)

	xlims = np.array(xlims)
	ylims = np.array(ylims) 

	# PLOTTING SOURCE POINTS
	if marker == "o": 	# 'o' used for open circle 
		ax.scatter(xsources, ysources, facecolor="none", edgecolor=color, s=size, label=label)
	elif marker == None: 	# default is a closed circle 
		ax.scatter(xsources, ysources, color=color, s=size, label=label)
	else: 
		ax.scatter(xsources, ysources, color=color, s=size, label=label, marker=marker)

	# PLOTTING POINT NAMES, IF GIVEN
	if labelpoints: 
		for i in range(len(labelpoints)): 
			ax.annotate(labelpoints[i], xy=[xsources[i], ysources[i]-.2], size=10, horizontalalignment="center")

	# Setting plot limits
	if not xlim:
		if invert_xlim: plt.xlim(max(xlims.T[0]), min(xlims.T[1]))
		else: plt.xlim(min(xlims.T[0]), max(xlims.T[1]))
	else: plt.xlim(xlim)
	if not ylim: 
		if invert_ylim: plt.ylim(max(ylims.T[0]), min(ylims.T[1]))
		else: plt.ylim(min(ylims.T[0]), max(ylims.T[1]))
	else: plt.ylim(ylim)

	# plotting mass track labels
	plt.xlabel(xlabel, labelpad=0, fontsize=fontsize)#, labelparams, labelpad=0)
	plt.ylabel(ylabel, labelpad=-10, fontsize=fontsize)#, labelparams, labelpad=-10)


	# If another title (such as name of the object) is given, plot
	# Adjusts size to make it fit the size of the plot
	if 0.6*figsize[0]*figsize[1] > 30: titlesize = 30
	else: titlesize = 0.6*figsize[0]*figsize[1]

	if not annotation_size: annotation_size = titlesize

	if title: ax.annotate(title, xy=(xlim[1]-0.35*abs(xlim[1]), ylim[1]+0.2*abs(ylim[1])), size=annotation_size)

	# If an annotation is given, add it to the bottom of the figure
	if annotation: ax.annotate(annotation, xy=(xlim[1]-0.05*abs(xlim[1]), ylim[0]-0.1*abs(ylim[1])), size=annotation_size, horizontalalignment="right")

	# if a subimage is given as a filename, read in. 
	if subimg and "." in subimg: 
		subimg = img.imread(subimg)

	# If an image is passed to overlay on plot, add 
	if hasattr(subimg, 'shape'): # tests if subimg was fed in
		try: 
			XY = [figsize[0], figsize[1]]
			ax2 = f.add_axes([.9 - 0.22*float(XY[1])/float(XY[0]), 0.66,  float(XY[1])/float(XY[0])*0.22, 0.22], zorder=1)
			ax2.imshow(subimg)
			ax2.axis('off')
		except: 
			im = plt.imread(subimg)
			XY = [figsize[0], figsize[1]]
			ax2 = f.add_axes([.9 - 0.22*float(XY[1])/float(XY[0]), 0.66,  float(XY[1])/float(XY[0])*0.22, 0.22], zorder=1)
			ax2.imshow(im)
			ax2.axis('off')

	# saving image, if prompted
	if savefile != None: save = True; pass;
	if save:
		if savefile == None: 
			savefile = df["ID"][0].split("X")[0] + "_" + xcolors + "_" + ycolors
		plt.savefig(savefile.split(".")[0]+".jpg", dpi=300, bbox_inches="tight")

	# Returning plot information, in case I need this later
	# need to retrieve ax if using both subimg followed by AddCMD

	cd(temp_dir)
	return f, ax
				

###-----------------------------------------------------------------------------------------------------

def AddCMD(df=None, xcolors=VI, ycolors=V, color="black", size=10, marker=None, label=None, f=None, ax=None, color_correction=[0,0]): 

	"""Adds multiple sets of points to a single CMD plot. Should be used after MakeCMD. If plots do not print as expected, call in f and ax from MakeCMD."""

	try: 
		xi = colors.index(xcolors)
		yi = colors.index(ycolors)


		# Making sure all of the necessary colors are available in the DataFrame
		try: 
			temp = (df[xcolors].values - color_correction[0]).tolist()
		except: 
			df_heads = df.columns.values.tolist()
			if xi == 3: df_xcol = (df[B].values - df[V].values - color_correction[0]).tolist()
			elif xi == 4: df_xcol = (df[V].values - df[I].values - color_correction[0]).tolist()
			elif xi == 5: df_xcol = (df[B].values - df[I].values - color_correction[0]).tolist()
			else: print(xcolors + " not available.")
			df_values = [df[i].values.tolist() for i in df_heads]
			df_heads.append(xcolors)
			df_values.append(df_xcol)
			df = BuildFrame(headers=df_heads, values=df_values)
			df_heads = df.columns.values.tolist()
		try: 
			temp = (df[ycolors] - color_correction[1]).values.tolist()
		except: 
			df_heads = (df.columns.values - color_correction[1]).tolist()
			if yi == 3: df_ycol = (df[B].values - df[V].values - color_correction[1]).tolist()
			elif yi == 4: df_ycol = (df[V].values - df[I].values - color_correction[1]).tolist()
			elif yi == 5: df_ycol = (df[B].values - df[I].values - color_correction[1]).tolist()
			else: print(ycolors + " not available.")
			df_values = [df[i].values.tolist() for i in df_heads]
			df_heads.append(ycolors)
			df_values.append(df_ycol)
			df = BuildFrame(headers=df_heads, values=df_values)
			df_heads = df.columns.values.tolist()


		if ax: 	# ax MUST be read in if subimg is used in MakeCMD. 
			if marker == "o": 	# 'o' used for open circle 
				ax.scatter(df[xcolors], df[ycolors], facecolor="none", edgecolor=color, s=size, label=label)
			elif marker == None: 	# default is a closed circle 
				ax.scatter(df[xcolors], df[ycolors], color=color, s=size, label=label)
		else: 
			if marker == "o": 	# 'o' used for open circle 
				plt.scatter(df[xcolors], df[ycolors], facecolor="none", edgecolor=color, s=size, label=label)
			elif marker == None: 	# default is a closed circle 
				plt.scatter(df[xcolors], df[ycolors], color=color, s=size, label=label)

	except: return "Failed to add points."		
	

###-----------------------------------------------------------------------------------------------------

def GetMasses(df, colors=[VI,V], setv="F547M", setb="F438W", seti="F814W"): 

	"""Assigns masses based on the position of the source on the CMD relative to the mass tracks. Outputs either a list of masses, or the DataFrame"""

	df = df.copy()

	# For calculating colors based on input xcolors, ycolors, and defined or default filters
	# used for mass tracks
	mass_colors = [setb, setv, seti]
	mass_colorpairs = [[setb, setv], [setv, seti], [setb, seti]]

	xi = colors.index(colors[0])   # the positions of the CMD colors lines
	yi = colors.index(colors[1])

	# Add a Mass column to the DataFrame, if one does not already exist
	try: 
		temp = [None]*len(df)
		df.insert(df.shape[1], "Mass", temp)
	except: pass;

	# Comparing each source to the mass tracks
	for i in range(len(df)): 
		for j in range(len(min_masses)): 
			track = min_masses[j]	# go through the tracks one by one
			# Looking for the correct filters to compare to
			try: 
				if xi > 2: xtemp = track[mass_colorpairs[xi-3][0]] - track[mass_colorpairs[xi-3][1]]
				else: xtemp = track[mass_colors[xi]]
				if yi > 2: track = track[mass_colorpairs[yi-3][0]] - track,[mass_colorpairs[yi-3][1]]
				else: ytemp = track[mass_colors[yi]]
			except: 
				xtemp = track[colors[xi]]
				ytemp = track[colors[yi]]
			# Interpolating the track so that we can compare magnitude at the xval of the source
			f1 = interp1d(xtemp, ytemp, kind="nearest", fill_value="extrapolate")
			linetemp = f1(df[colors[0]][i])
			if df[colors[1]][i] < linetemp: 
				df["Mass"][i] = mass_labels[j]
				pass;
			elif df[colors[1]][i] == linetemp: 
				df["Mass"][i] = mass_labels[j+1]
				break;
			else:  
				df["Mass"][i] = mass_labels[j]
				break;

	return df


###-----------------------------------------------------------------------------------------------------

def GetPhots(filename=None, ids=None, outfile=None, savephot=None, getaps=False, getcoords=False, aperture=None): 

	"""Retrieves the photometry from a given .ph file. If IDs are given, all other objects are filtered out. Outfile saves the resulting photometry."""

	if filename == None: filename = input("Photometry file name?: ")
	phot = np.genfromtxt(filename, dtype=float)

	# If IDs are given as a file or a list, filter out all other phots
	if ids: 
		try: 
			ids = np.genfromtxt(ids, dtype=float)
			ids = [int(i) for i in ids]
		except: pass;
		phot = np.array([phot[int(i)] for i in ids])		

	if savephot: 
		f = open(savephot, "wb")
		np.savetxt(f, phot)
		print(savephot + " saved!")

	# Getting the number of apertures
	numaps = int((phot.shape[-1] - 2)/3)
	
	if getaps: 
		if len(phot.shape) > 1: aps = phot.T[2:numaps+2].T[0]
		else: aps = phot[2:numaps+2]
	if getcoords: 
		coords = []
		for i in range(len(phot)): 
			coords.append([phot.T[0][i], phot.T[1][i]])

	phot = phot.T[2+numaps:-1*numaps]
	phot = phot.T

	if outfile: 
		f = open(outfile, "wb")
		np.savetxt(f, phot)
		print(outfile + " saved!")
	
	if aperture: phot = np.array(phot).T[default_aps.index(aperture)].tolist()

	if getaps and getcoords: 
		return phot, coords, aps
	elif getaps: 
		return phot, aps
	elif getcoords: 
		return phot, coords
	else: return phot


###-----------------------------------------------------------------------------------------------------

def CorrectAp(filename=None, goodstars=None, verbose=True, getstd=False, savefile=None, savecorr=None, apertures=[3,20], savecoords=None, colors=None, labels=None, randomize=False, nrand=100, saverand=None, overwrite_rand=True, field=None, xlim=[0,20], ylim=[26,18]): 

	"""For calculating the aperture correction using a sample of previously chosen and measured, bright, isolated stars. Parameter 'goodstars' may be a list of good stars ids to pull from the file listed in 'filename' or a .txt file with the ids listed. """
	
	# savecoords can be bool or list. If given, getcoords=True
	if savecoords: getcoords = True
	else: savecoords=False; getcoords = False

	# If filename is given as a list, treat these as a set of files to compare
	# Needs to convert filename as a list if only one file given. 
	if not isinstance(filename, list): filename = [filename]
	phots = [] # contains photometry of all stars in .ph file(s).
	aps = []   # list of aperture sizes from the file
	if getcoords: coords = []	# might not need coords
	
	# If goodstars is given, either pull list from file or read in for GetPhots
	# This is overwritten if randomize is set to True
	if goodstars: 
		if ".txt" in goodstars: goodstars = np.genfromtxt(goodstars,dtype=int).tolist()

	# NOTE: This is written as if more than one filename can be passed, but that doesn't work well in principle. 
	#       Due to this, goodstars is assumed to be a single file, not a list of files, 
	#       and filename should NOT be read in as a list (not fixed as of April 7, 2021)
	for i in filename:
		try: 
			temp_phots, temp_coords, temp_aps = GetPhots(filename=i, ids=goodstars, getaps=True, getcoords=getcoords)
			phots.append(temp_phots)
			coords.append(temp_coords)
			aps.append(temp_aps)
		except: 
			temp_phots, temp_aps = GetPhots(filename=i, ids=goodstars, getaps=True)
			phots.append(temp_phots)
			aps.append(temp_aps)

		#--- END FOR LOOP


	# If randomize = True, randomly draw from the input .ph file to test the correction on. 
	# This is helpful for when finding the perfect stars is difficult. 
	# If a good sample of perfect stars isn't picked up by this method, run again. 
	if randomize: 
		for i in range(len(phots)): # for each .ph file contained in phots...
			temp = phots[i]     # look at the photometry in the current field
			# Select nrand number of random stars
			temp_select_temp = random.sample(range(0, len(temp)), nrand)
			temp_select_temp.sort()
			temp_select = []

			cont = True # counter for continuing while loop

			# Checks the validity of the radial profile.
			# Stars that dip or don't flatten should be automatically excluded.
			# While the length of temp_select is less than nrand, add more stars.
			while cont: 
				if len(temp_select) < nrand:
					for j in temp_select_temp: 
						temp2 = temp[j] # full radial profile of chosen random star
						# checking for dips and flattened ends
						if temp2[2] < 25.5 and temp2[-1] > 19:
							if all(k>l for k,l in zip(temp2, temp2[1:])): # ensures the profile always decreases
								if all(k-l<0.15 for k,l in zip(temp2[3:], temp2[4:])):
									if j not in temp_select: temp_select.append(j)
					temp_select_temp = random.sample(range(0, len(temp)), nrand-len(temp_select))
				else: cont = False; pass;

			temp_select.sort()
			phots[i] = phots[i][temp_select]

	# Get the magnitudes and apertures from the file
	# Plot all of the stars magnitudes vs. aperture radius to compare straightness
	# Plot each star mag vs. ap *** (start of loop)
	# Ask for which star is bad
	# Invert to create mask for "good stars"
	# Plot all "good stars"
	# As if process needs to be repeated *** (reset loop)
	# If so, start from plotting each star and repeat until good
	# Take the median and std of the aperture corrections between given apertures (?)
	# Added 'compare=None' to give me the option to compare same stars in different filters
	#	and remove star across the board (updated: 7/22/19)

	# Getting apertures to compare
	compare = []
	for i in aps: compare.append([i.tolist().index(apertures[0]),i.tolist().index(apertures[1])])  
	print("Comparing apertures: " + str(int(aps[0][compare[0][0]])) + " and " + str(int(aps[0][compare[0][1]])))

	# Setting up files for saving coordinates
	if savecoords and isinstance(savecoords, bool): 
		# if savecoords indicated but not given, use the filenames as savecoords names
		savecoords = []
		for i in filename: 
			savecoords.append(re.split(".ph", i)[0] + ".coords")
	elif savecoords and not isinstance(savecoords, list): 
		savecoords = [savecoords]

	if savefile and not isinstance(savefile, list): 
		savefile = [savefile]
	

	# Setting up colors and labels for plotting, if needed
	if colors and isinstance(colors, list): 
		while len(colors) != len(filename): 
			print("Number of colors does not match number of files (" + str(len(filename)) + " needed)") 
			colors = input("Plot line colors (separate by commas): ")
			colors = re.split("\W+", colors)
	if labels and isinstance(labels, list): 
		while len(labels) != len(filename): 
			print("Number of labels does not match number of files (" + str(len(filename)) + " needed)") 
			labels = input("Plot line lables (separate by commas): ")
			labels = re.split(",", labels)	

	# Convert everything to arrays for easy manipulation
	phots = np.array(phots)
	aps = np.array(aps)
	if getcoords: coords = np.array(coords)

	# Plotting all stars to see if all crooked ones are removed
	for i in range(len(phots)):
		for j in phots[i]:  
			plt.plot(aps[i], j)
			plt.ylim(ylim[0],ylim[1])
			plt.xlim(xlim[0],xlim[1])
		plt.title(filename[i])
		plt.show()

	try: 
		userin = input("Remove stars?: ").lower()
		if userin[0] == "y" or userin[0] == "r": 
			repeat = True
		else: repeat = False; pass;
	except: repeat = False; pass;

	while repeat == True: 
		# If verbose is set to true, then ask if star should be removed after each plot. 
		if verbose == True: 

			bad = []	# keeping track of bad stars

			# Plotting all given photometry per single star with object number
			for j in range(len(phots[0])): # for each star in files
				for i in range(len(phots)):	# for all files
					if colors and labels: 
						plt.plot(aps[i], phots[i][j], color=colors[i], label=labels[i])
						plt.legend()
					elif colors: plt.plot(aps[i], phots[i][j], color=colors[i]) 
					elif labels: 
						plt.plot(aps[i], phots[i][j], label=labels[i]) 
						plt.legend()
					else: plt.plot(aps[i], phots[i][j])
					plt.ylim(ylim)
					plt.xlim(xlim)
					if randomize: plt.ylabel("Star ID: " + str(temp_select[j]))  # So that I can keep track of the 'good' stars

				plt.title("Star " + str(j), size=20)
				plt.show()	# plot all files per star on one plot

				# After plotting, ask if it should be removed (as per verbose)
				# This just removes star # j
				try: 
					userin = input("Remove?: ").lower()
					if userin[0] == "y" or userin[0] == "r": bad.append(j)
					else: plt.close(); pass;
				except: plt.close(); pass;

			#--- END IF STATEMENT

		# If verbose is set to false, plot each star first, then ask for list of stars to remove.
		else: 
			# Plotting each star with object number
			for j in range(len(phots[0])): 
				for i in range(len(phots)):
					if colors and labels: plt.plot(aps[i], phots[i][j], color=colors[i], label=labels[i])
					elif colors: plt.plot(aps[i], phots[i][j], color=colors[i]) 
					elif labels: plt.plot(aps[i], phots[i][j], label=labels[i]) 
					else: plt.plot(aps[i], phots[i][j])
					plt.ylim(ylim[0],ylim[1])
					plt.xlim(xlim[0],xlim[1])
				plt.title(j, size=20)
				plt.show()

			# Removing bad stars from the list
			bad = input("Remove which star number? (separate by commas): ")

			if len(bad) > 0: 
				bad = re.split("\W+", bad)
				bad = [int(x) for x in bad]

			#--- END ELSE STATEMENT

		print("\nRemoving bad stars...\n")
		
		temp_phots =[]

		for i in phots: # looking at all stars in each file
			temp_phots.append([x for j,x in enumerate(i) if not j in bad])
		phots = np.array(temp_phots)

		if getcoords: # do the same if coords are given
			temp_coords = []
			for i in coords: 
				temp_coords.append([x for j,x in enumerate(i) if not j in bad])
			coords = np.array(temp_coords)
		else: pass;

		# And if randomize is set, remove bad stars from the randomized list.
		if randomize: 
			temp_select = [temp_select[i] for i in range(len(temp_select)) if not i in bad]
		else: pass;


		# Plotting all stars to see if all crooked ones are removed
		for i in range(len(phots)):
			for j in phots[i]:  
				plt.plot(aps[i], j)
				plt.ylim(ylim[0],ylim[1])
				plt.xlim(xlim[0],xlim[1])
			plt.title(filename[i])
			plt.show()

		# Testing current aperture correction and stds
		print("Aperture corrections and (standard deviations)") 
		for i in range(len(phots)): 
			print(str(nanmedian(np.array(phots)[i].T[compare[i][0]] - np.array(phots)[i].T[compare[i][1]])) + " ("  + str(nanstd(np.array(phots)[i].T[compare[i][0]] - np.array(phots)[i].T[compare[i][1]])) + ")")

		try: 	# remove more if currentaperture correction and std isn't good enough
			userin = input("Remove more?: ").lower()
			if userin[0] == "y": repeat = True; pass;
			elif userin == "reset" or userin == "restart": 	# reset removal if needed
				repeat = True
				print("\nResetting aperture corrections...\n")
				phots = []
				aps = []
				if getcoords: coords = []	# might not need coords
				for i in filename:
					try: 
						temp_phots, temp_coords, temp_aps = GetPhots(filename=i, getaps=True, getcoords=getcoords)
						phots.append(temp_phots)
						coords.append(temp_coords)
						aps.append(temp_aps)
					except: 
						temp_phots, temp_aps = GetPhots(filename=i, getaps=True)
						phots.append(temp_phots)
						aps.append(temp_aps)

				# Printing all photometries per file
				for i in range(len(phots)):
					for j in phots[i]:  
						plt.plot(aps[i], j)
						plt.ylim(ylim[0],ylim[1])
						plt.xlim(xlim[0],xlim[1])
					plt.title(filename[i])
					plt.show()
			else: repeat = False; break;
		except: repeat = False; break; # If no other bad stars exist, break the while cycle
		
		# Converting data into arrays for easier handling
		phots = np.array(phots)
		try: coords = np.array(coords)
		except: pass;
	
		#--- END WHILE LOOP

	# Saving results to a file if filename is given
	if savefile:
		for i in range(len(savefile)): 
			try:
				f = open(savefile[i], "wb")
				np.savetxt(f, phots[i]) 
				f.close()
			except: print("\nSaving " + savefile[i] + " failed.")
	if savecoords: 
		for i in range(len(savecoords)):
			try:
				f = open(savecoords[i], "wb")
				np.savetxt(f, np.column_stack([coords[i].T[0], coords[i].T[1]]))
				f.close()
			except: print("\nSaving " + savecoords[i] + " failed.")

	# Currently only allows for a single random star selection save
	# If the file exists, add new stars to bottom. If it doesn't, create it.

	if saverand: 
		if not overwrite_rand:
			f = open(saverand, "a+")
			np.savetxt(f, np.column_stack([np.array(temp_select).T]))
			f.close()
			print("Added " + str(len(temp_select)) + " good stars to " + saverand)
		else: 
			f = open(saverand, "w+")
			np.savetxt(f, np.column_stack([np.array(temp_select).T]))
			f.close()
			print("Saved good stars as " + saverand)

	print("Done")

	# Returning results
	return_phots = []
	for i in range(len(phots)): 
		return_phots.append(nanmedian(phots[i].T[compare[i][0]] - phots[i].T[compare[i][1]]))
	
	return_std = []
	for i in range(len(phots)): 
		return_std.append(nanstd(phots[i].T[compare[i][0]] - phots[i].T[compare[i][1]]))
	
	if savecorr: 
		f = open(savecorr, "a+")
		f.write(str(field)+" "+str(return_phots[0])+" "+str(return_std[0])+"\n")
		f.close()

	if getstd: return return_phots, return_std
	else: return return_phots


###-----------------------------------------------------------------------------------------------------

def CorrectMags(frame=None, phots=None, corrections=None, field=None, apertures=[3,20], headers=[V, B, I], instrument="ACS", filters=["F606W", "F435W", "F814W"], distance=3.63e6, savefile=None, ID_header=ID, coord_headers=[X, Y], extinction=[0,0,0,0]): 

	"""Calculating magnitudes with given aperture corrections. The input 'instrument' can be 'ACS' or 'WFC3', which defines which EEF file to read from. Filters should be read in the order [V,B,I]. If given, 'extinction' should also be in [Av,Ab,Ai,Au] order. (NOTE: note RGB) or [V,B,I,U], if U is given. Corrections should also be read in VBI order. """

	try: frame = frame.copy()
	except: pass;

	# If U is given, add U corrections to all commands
	if len(filters) == 4: U_true = True
	else: U_true = False
	
	current_dir = pwd()
	# Loading in the EEFs file
	cd("/home/qiana/Documents/Research/XRB/")
	EEFs = LoadSources(instrument + "_EEFs.txt", verbose=False)
	cd(current_dir)

	# Reading in the 20px (or highest given) EEF from each filter 
	V_EEF = Find(EEFs, "Filter = " + filters[0])[str(apertures[-1])][0]
	B_EEF = Find(EEFs, "Filter = " + filters[1])[str(apertures[-1])][0]
	I_EEF = Find(EEFs, "Filter = " + filters[2])[str(apertures[-1])][0]
	if U_true: U_EEF = Find(EEFs, "Filter = " + filters[3])[str(apertures[-1])][0]

	#EEFracts=[Vee["20"],Bee["20"],Iee["20"],Uee["20"]]
	#V_EEF = EEFracts[0]
	#B_EEF = EEFracts[1]
	#I_EEF = EEFracts[2]
	#U_EEF = EEFracts[3]

	if not phots: 
		phots = []
		phots.append(frame[headers[0]])
		phots.append(frame[headers[1]])
		phots.append(frame[headers[2]])
		if U_true: phots.append(frame[headers[3]])

	dmod = 5.*log10(distance/10.)

	# Converting to absolute magnitudes? 
	V_M = phots[0] - dmod - extinction[0]
	B_M = phots[1] - dmod - extinction[1]
	I_M = phots[2] - dmod - extinction[2]
	if U_true: U_M = phots[3] - dmod - extinction[3]

	# Finding the proper corrections factor. If none given, find proper defaults.
	if not corrections: 
		if not field: 
			try: field = frame["Field"].values.tolist()
			except: field = input("Image Field? (f_):")
		if isinstance(field, list): 
			corrections = [[Corrections[V][int(i)-1], Corrections[B][int(i)-1], Corrections[I][int(i)-1], Corrections[U][int(i)-1]] for i in field]
			corrections = np.array(corrections)
		else: 
			try: corr = int(re.split("f", field.lower())[-1][0])
			except: corr = int(field)
			corrections = [Corrections[V][corr-1], Corrections[B][corr-1], Corrections[I][corr-1], Corrections[U][corr-1]]
			corrections = np.array(corrections)

	# Calculating the corrections
	try: 
		V_corr = np.array(V_M) - corrections.T[0] + 2.5*log10(float(V_EEF))
		B_corr = np.array(B_M) - corrections.T[1] + 2.5*log10(float(B_EEF))
		I_corr = np.array(I_M) - corrections.T[2] + 2.5*log10(float(I_EEF))
		if U_true: U_corr = np.array(U_M) - corrections.T[3] + 2.5*log10(float(U_EEF))
	except: 
		V_corr = np.array(V_M) - corrections[0] + 2.5*log10(float(V_EEF))
		B_corr = np.array(B_M) - corrections[1] + 2.5*log10(float(B_EEF))
		I_corr = np.array(I_M) - corrections[2] + 2.5*log10(float(I_EEF))
		if U_true: U_corr = np.array(U_M) - corrections[3] + 2.5*log10(float(U_EEF))

	for i in range(len(V_corr)): 
		if V_corr[i] > 100 or V_corr[i] < -100: V_corr[i] = np.nan
		if B_corr[i] > 100 or B_corr[i] < -100: B_corr[i] = np.nan
		if I_corr[i] > 100 or I_corr[i] < -100: I_corr[i] = np.nan
		if U_true: 
			if U_corr[i] > 100 or U_corr[i] < -100: U_corr[i] = np.nan

	if savefile: 
		try: Mags = BuildFrame(headers=[ID_header, X, Y, V, B, I, U, VI, BV, BI], \
		         values = [frame[ID_header], frame[coord_headers[0]], frame[coord_headers[1]], V_corr, B_corr, I_corr, U_corr, V_corr-I_corr, B_corr-V_corr, B_corr-I_corr])
		except: 
			Mags = BuildFrame(headers=[ID_header, X, Y, V, B, I, U, VI, BV, BI], \
		         values = [frame[ID_header], frame[coord_headers[0]], frame[coord_headers[1]], V_corr, B_corr, I_corr, V_corr-I_corr, B_corr-V_corr, B_corr-I_corr])

		SaveDF(Mags, savefile)

	try: 
		frame[headers[0]] = V_corr
		frame[headers[1]] = B_corr
		frame[headers[2]] = I_corr
		if U_true: frame[headers[3]] = U_corr
		return frame
	except: 
		try: return V_corr, B_corr, I_corr, U_corr
		except: return V_corr, B_corr, I_corr

###-----------------------------------------------------------------------------------------------------

###-----------------------------------------------------------------------------------------------------

def CorrectMag(frame=None, phots=None, correction=None, field=None, apertures=[3,20], instrument="ACS", filt="F606W", distance=3.63e6, savefile=None, ID_header=ID, coord_headers=[X, Y], extinction=0): 

	"""Calculating magnitude with given aperture correction, like CorrectMags, but specifically for a single input filter (so that it doesn't require all filters to be given if only one measurement is needed). The input 'instrument' can be 'ACS' or 'WFC3', which defines which EEF file to read from. Filters should be read in the order [V,B,I]. If given, 'extinction' should also be in [Av,Ab,Ai,Au] order. (NOTE: note RGB) or [V,B,I,U], if U is given. Corrections should also be read in VBI order. """

	try: frame = frame.copy()
	except: pass;

	current_dir = pwd()
	# Loading in the EEFs file
	cd("/home/qiana/Documents/Research/XRB/")
	EEFs = LoadSources(instrument + "_EEFs.txt", verbose=False)
	cd(current_dir)

	if filt in V_filts: 
		EEF = Find(EEFs, "Filter = " + filt)[str(apertures[-1])][0]
		header = V
	elif filt in B_filts: 
		EEF = Find(EEFs, "Filter = " + filt)[str(apertures[-1])][0]
		header = B
	elif filt in I_filts: 
		EEF = Find(EEFs, "Filter = " + filt)[str(apertures[-1])][0]
		header = I
	elif filt in U_filts: 
		EEF = Find(EEFs, "Filter = " + filt)[str(apertures[-1])][0]
		header = U

	if not phots: 
		phots = []
		phots.append(frame[header])

	dmod = 5.*log10(distance/10.)

	# Converting to absolute magnitudes? 
	Mag = phots - dmod - extinction

	# Finding the proper corrections factor. If none given, find proper defaults.
	# Double check this. I believe this is now obsolete
	#if not correction: 
	#	if not field: 
	#		try: field = frame["Field"].values.tolist()
	#		except: field = input("Image Field? (f_):")
	#	if isinstance(field, list): 
	#		correction = [Corrections[header][int(i)-1] for i in field]
	#		correction = np.array(correction)
	#	else: 
	#		try: corr = int(re.split("f", field.lower())[-1][0])
	#		except: corr = int(field)
	#		correction = [Corrections[header][corr-1]]
	#		correction = np.array(correction)

	# Calculating the corrections
	try: 
		corr = np.array(Mag) - correction + 2.5*log10(float(EEF))
	except: 
		corr = np.array(Mag) - correction + 2.5*log10(float(EEF))

	#for i in range(len(corr)): 
	#	if corr[i] > 100 or corr[i] < -100: corr[i] = np.nan

	if savefile: 
		try: Mags = BuildFrame(headers=[ID_header, X, Y, header], \
		         values = [frame[ID_header], frame[coord_headers[0]], frame[coord_headers[1]], corr])
		except: 
			Mags = BuildFrame(headers=[ID_header, X, Y, header], \
		         values = [frame[ID_header], frame[coord_headers[0]], frame[coord_headers[1]], corr])

		SaveDF(Mags, savefile)

	try: 
		frame[header] = corr.tolist()[0]
		return frame
	except: 
		return corr

###-----------------------------------------------------------------------------------------------------


def CalColors(phots=None, reddening=[0,0,0]):	
	
	"""Calculates the colors from the given VBI magnitudes. Reddening may be used to add a reddening correction in the order [E(B-V), E(V-I), E(B-I)], and extinction should be given as [V,B,I]. NOTE: still need to add the possibility of U."""
	
	if not isinstance(phots, list): 
		phots = [phots[V].values.tolist(), phots[B].values.tolist(), phots[I].values.tolist()]


	try: 
		vi = (np.array(phots[0] - phots[2]) - reddening[1]).tolist()
		bv = (np.array(phots[1] - phots[0]) - reddening[0]).tolist()
		bi = (np.array(phots[1] - phots[2]) - reddening[2]).tolist()
	except: 
		vi = [phots[0][i] - phots[2][i] - reddening[1] for i in range(len(phots[0]))]
		bv = [phots[1][i] - phots[0][i] - reddening[0] for i in range(len(phots[0]))]
		bi = [phots[1][i] - phots[2][i] - reddening[2] for i in range(len(phots[0]))]


	print("Returning VI, BV, and BI.")
	
	return vi, bv, bi


###-----------------------------------------------------------------------------------------------------

def CorrectApOld(filename=None, verbose=True, getstd=False, savefile=None, apertures=[3,20], savecoords=None, compare=None): 

	"""For calculating the aperture correction."""
	
	if savecoords: getcoords = True
	else: getcoords = False

	try: phots, coords, aps = GetPhots(filename=filename, getaps=True, getcoords=getcoords)
	except: phots, aps = GetPhots(filename=filename, getaps=True)

	# Get the magnitudes and apertures from the file
	# Plot all of the stars magnitudes vs. aperture radius to compare straightness
	# Plot each star mag vs. ap *** (start of loop)
	# Ask for which star is bad
	# Invert to create mask for "good stars"
	# Plot all "good stars"
	# As if process needs to be repeated *** (reset loop)
	# If so, start from plotting each star and repeat until good
	# Take the median and std of the aperture corrections between given apertures (?)
	# Added 'compare=None' to give me the option to compare same stars in different filters
	#	and remove star across the board (updated: 7/22/19)

	# Getting apertures to compare
	compare = [aps.tolist().index(apertures[0]),aps.tolist().index(apertures[1])]  
	print("Comparing apertures: " + str(int(aps[compare[0]])) + " " + str(int(aps[compare[1]]))) 

	# Plotting all stars to see if all crooked ones are removed
	plt.figure(figsize=(8,6))
	for i in phots: 
		plt.plot(aps, i)
		plt.ylim(26, 18)
	plt.show()

	try: 
		if input("Remove stars?: ").lower()[0] == "y" or input("Remove stars?: ").lower()[0] == "r": 
			repeat = True
		elif se: repeat = False; pass;
	except: repeat = False; pass;

	# If verbose is set to true, then ask if star should be removed after each plot. 
	# UNDER CONSTRUCTION 
	if verbose == True: 
		while repeat == True: 
			bad = []

			# Plotting each star with object number
			for i in range(len(phots)): 
				j = phots[i]
				plt.plot(aps, j)
				plt.title(i, size=20)
				plt.ylim(26, 18)
				plt.show()

				# After plotting, ask if it should be removed (as per verbose)
				try: 
					if input("Remove?: ").lower()[0] == "y" or input("Remove?: ").lower()[0] == "r": bad.append(i)
					else: plt.close(); pass;
				except: plt.close(); pass;

			print("\nRemoving bad stars...\n")
			phots = [x for i,x in enumerate(phots) if not i in bad]
			try: coords = [x for i,x in enumerate(coords) if not i in bad]
			except: pass;

			# Plotting all stars to see if all crooked ones are removed
			plt.figure(figsize=(10,8))
			for i in phots: 
				plt.plot(aps, i)
				plt.ylim(26, 18)
			plt.show()
			print("Aperture correction: " + str(nanmedian(np.array(phots).T[compare[0]] - np.array(phots).T[compare[1]])) + "\nStandard deviation: "  + str(nanstd(np.array(phots).T[compare[0]] - np.array(phots).T[compare[1]])))
			try: 
				temp = input("Remove more?: ").lower()
				if temp[0] == "y": repeat = True; pass;
				elif temp == "reset" or temp == "restart": 
					repeat = True
					print("\nResetting aperture corrections...\n")
					try:  
						phots, coords, aps = GetPhots(filename=filename, getaps=True, getcoords=getcoords)
					except: phots, aps = GetPhots(filename=filename, getaps=True)
					for i in phots: 
						plt.plot(aps, i)
						plt.ylim(26, 18)
					plt.show()
			except: repeat = False; break; # If no other bad stars exist, break the while cycle

	else: 
		# If verbose is set to false, plot each star first, then ask for list of stars to remove.
		# Repeat plotting magnitudes vs. apertures until all bad stars are removed
		while repeat == True: 

			# Plotting each star with object number
			for i in range(len(phots)): 
				j = phots[i]
				plt.plot(aps, j)
				plt.title(i, size=20)
				plt.ylim(26, 18)
				plt.show()

			# Removing bad stars from the list
			bad = input("Remove which star number? (separate by commas): ")

			if len(bad) > 0: 
				bad = re.split("\W+", bad)
				bad = [int(x) for x in bad]
				phots = [x for i,x in enumerate(phots) if not i in bad]
				try: coords = [x for i,x in enumerate(coords) if not i in bad]
				except: pass;

				# Plotting all stars to see if all crooked ones are removed
				plt.figure(figsize=(10,8))
				for i in phots: 
					plt.plot(aps, i)
					plt.ylim(26, 18)
				plt.show()
				print("Aperture correction: " + str(nanmedian(np.array(phots).T[compare[0]] - np.array(phots).T[compare[1]])) + "\nStandard deviation: "  + str(nanstd(np.array(phots).T[compare[0]] - np.array(phots).T[compare[1]])))
				try: 
					temp = input("Remove more?: ").lower()
					if temp[0] == "y": repeat = True
					elif temp == "reset" or temp == "restart": 
						repeat = True
						print("\nResetting aperture corrections...\n")
						try: 	
							phots, coords, aps = GetPhots(filename=filename, getaps=True, getcoords=getcoords)
						except:  
							phots, aps = GetPhots(filename=filename, getaps=True)
							plt.figure(figsize=(10,8))
						for i in phots: 
							plt.plot(aps, i)
							plt.ylim(26, 18)
						plt.show()
				except: repeat = False; pass;
			else: repeat = False; break; # If no other bad stars exist, break

	# Converting data into arrays for easier handling
	phots = np.array(phots)
	coords = np.array(coords)

	# Saving results to a file if filename is given
	if savefile:
		for i in savefile: 
			try:
				f = open(savefile[i], "wb")
				np.savetxt(f, phots[i]) 
				f.close()
			except: print("\nSaving " + savefile[i] + " failed.")
	if savecoords: 
		for i in savecoords:
			try:
				f = open(savecoords[i], "wb")
				np.savetxt(f, np.column_stack([coords[i].T[0], coords[i].T[1]]))
				f.close()
			except: print("\nSaving " + savecoords[i] + " failed.")

	# Returning results
	if getstd == False: # Check if returning std
		return nanmedian(phots.T[compare[0]] - phots.T[compare[1]])
	else: 
		return nanmedian(phots.T[compare[0]] - phots.T[compare[1]]), nanstd(phots.T[compare[0]] - phots.T[compare[1]])

	

###-----------------------------------------------------------------------------------------------------

def GetTracks(): 

	"""Reads in info for tracks from the .txt or .dat files output from Padova Mass Tracks website."""



				

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

# Pull directory where XRBID files are saved
file_dir = os.path.dirname(os.path.abspath(__file__))
curr_dir = pwd()

from XRBID.DataFrameMod import Find, BuildFrame
from XRBID.Sources import LoadSources

default_aps = [i for i in range(1,31)] #[0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.]

# Labels associated with each line. Sources BELOW each line (but above the previous) get this label
mass_labels = ["Low", "Intermediate", "Intermediate", "High"]

# Old code refers to headers as defined below. These will need to be updated - 4/18/25
V = "V"
B = "B" 
I = "I"
U = "U"

# Calling in tracks for future use by MakeCMD
cd(file_dir)

wfc3_masses = [pd.read_csv("isoWFC3_1Msun.frame"), pd.read_csv("isoWFC3_3Msun.frame"), pd.read_csv("isoWFC3_5Msun.frame"),
	       pd.read_csv("isoWFC3_8Msun.frame"), pd.read_csv("isoWFC3_20Msun.frame")]

acs_masses = [pd.read_csv("isoACS_WFC_1Msun.frame"), pd.read_csv("isoACS_WFC_3Msun.frame"), pd.read_csv("isoACS_WFC_5Msun.frame"),
	      pd.read_csv("isoACS_WFC_8Msun.frame"), pd.read_csv("isoACS_WFC_20Msun.frame")]

cd(curr_dir)

###-----------------------------------------------------------------------------------------------------

def MakeCMD(sources=False, xcolor=None, ycolor=None, xmodel=None, ymodel=None, figsize=(6,4), xlim=None, ylim=None, color="black", size=10, marker=None, label=None, save=False, savefile=None, title=None, subimg=None, annotation=None, annotation_size=None, imshow=True, fontsize=15, shift_labels=[[0,0],[0,0],[0,0],[0,0],[0,0]], set_labels=None, instrument="ACS", color_correction=[0,0], labelpoints=False, file_dir=False): 

	"""Makes a CMD from a given set of points, either from a list or an input dataframe.

	PARAMETERS: 
	sources 	[pd.dataframe, list]: 	Input may either be a pandas dataframe containing the appropriate magnitudes required for the CMD, 
				       		or a list of coordinates in the format [[xs],[ys]]. 
	xcolor 		[str or list]: 		The name of the color or magnitude to plot in the x-axis, as it is labeled in the dataframe. 
						This will be the default x-axis label. If a list is given, it is assumed to contain the 
						column names of the input dataframe to be subtracted (e.g. F555W - F814W)
	ycolor 		[str or list]: 		The name of the color or magnitude to plot in the y-axis, as it is labeled in the dataframe. 
						This will be the default y-axis label. If a list is given, it is assumed to contain the column 
						names of the input dataframe to be subtracted (e.g. F555W - F814W)
	xmodel 		[str or list]: 		The magnitude(s) of the filter(s) to be used from the stellar models for the x-axis. 
						If given as a list, it is assumed the color is xmodel[0] - xmodel[1] 
	ymodel 		[str or list]: 		The magnitude(s) of the filter(s) to be used from the stellar models for the y-axis. 
						If given as a list, it is assumed the color is ymodel[0] - ymodel[1] 
	figsize 	[tuple] (6,4): 		The desired dimensions of the figure. 
	xlim 		[tuple] (None):		The limits on the x-axis. If none are given, the limits are assumed to be the (xmin - 1, xmax + 1) 
	ylim 		[tuple] (None):		The limits on the y-axis. If none are given, the limits are assumed to be (ymax + 1, ymin - 1)
	color 		[str] ("black"): 	Marker color
	size 		[int] (10): 		Marker size
	marker 		[str] (None):		The style of the marker. Defaults to a filled point. If "o" is given, 
						the marker is set to be an open circle. 
	label 		[str] (None):		The legend label to assign to the input points from sources. 
	save 		[bool] (False): 	Sets whether to automatically same the CMD image. 
	savefile 	[str] (None): 		The name to assigned the saved CMD image. 
	title 		[str)] (None): 		Title of the figure, to be placed near, but not at, the top of the figure.
	subimg 		[str] (None): 		Filename of an image to include in the corner of the CMD. 
						This is to allow a subplot of the XRB plotted to be shown within the CMD. 
	annotation 	[str] (None): 		Additional annotation to add to the bottom corner of the CMD (usually XRB ID)
	annotation_size [int] (None): 		Annotation fontsize 
	imshow 		[bool] (True): 		Shows the plot.
	fontsize 	[int] (20): 		The fontsize of text other than the annoations
	shift_labels 	[list]: 		List of x and y coordinate distances by which to shift each of the model mass labels.
						Defaults to [[0,0],[0,0],[0,0],[0,0],[0,0]] (no shifts)
	set_labels 	[list] (None):		Sets the position of labels. If none is given, positions are automatically calculated.
	instrument 	[str] ("ACS"): 		Name of the instrument used, to determine which models to call. 
	color_correction [list] ([0,0]):	Corrections on the x and y position of the sources. Defaults to no correction. 
	labelpoints 	[list] (False): 	Labels to add to each point. If none are given, defaults to False and no labels added. 
	file_dir 	[str]: 			The directory within which the models may be found. By default, the 
						code attempts to find this automatically, but if it fails, it will 
						prompt the user to input the directory manually. 

	RETURNS: 
	f, ax: 		Arguments defining the figure, which can be used to add more points to the CMD after the initial plotting.
 
"""

	# Setting the style of my plots
	#fontparams = {'font.family':'stix'}
	#labelparams = {'family':'stix', 'size':fontsize}
	
	curr_dir = pwd()

	# If no file directory is given, assume the files we need are in the same directory
	# where the module is saved
	#if not file_dir: 
	#	file_dir = os.path.dirname(os.path.abspath(__file__))

	#try: cd(file_dir)
	#except: print("Directory containg CMD models not found.\nPlease check and input the correct directory manually with file_dir.")

	# Reading in the appropriate models based on the instrument given.
	#if instrument.upper() =="WFC3":
	#	mass1 = pd.read_csv("isoWFC3_1Msun.frame")
	#	mass3 = pd.read_csv("isoWFC3_3Msun.frame")
	#	mass5 = pd.read_csv("isoWFC3_5Msun.frame")
	#	mass8 = pd.read_csv("isoWFC3_8Msun.frame")
	#	mass20 = pd.read_csv("isoWFC3_20Msun.frame")
	#elif instrument.upper() =="ACS": 
	#	mass1 = pd.read_csv("isoACS_WFC_1Msun.frame")
	#	mass3 = pd.read_csv("isoACS_WFC_3Msun.frame")
	#	mass5 = pd.read_csv("isoACS_WFC_5Msun.frame")
	#	mass8 = pd.read_csv("isoACS_WFC_8Msun.frame")
	#	mass20 = pd.read_csv("isoACS_WFC_20Msun.frame")

	if instrument.upper() =="WFC3": masses = wfc3_masses # list of DataFrames of each mass model
	else: masses = acs_masses
	mass_labels = ["1 M$_\odot$", "3 M$_\odot$", "5 M$_\odot$", "8 M$_\odot$", "20 M$_\odot$"]

	#cd(curr_dir) 

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
	if sources: 
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
	if sources: 
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

	cd(curr_dir)
	return f, ax
				

###-----------------------------------------------------------------------------------------------------

def AddCMD(df=None, xcolor=False, ycolor=False, color="black", size=10, marker=None, label=None, f=None, ax=None, color_correction=[0,0]): 

	"""Adds multiple sets of points to a single CMD plot. Should be used after MakeCMD. 
	If plots do not print as expected, call in f and ax from MakeCMD.
	NOTE: This code us currently under construction"""

	
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

	try: 
		if ax: 	# ax MUST be read in if subimg is used in MakeCMD. 
			if marker == "o": 	# 'o' used for open circle 
				ax.scatter(xsources, ysources, facecolor="none", edgecolor=color, s=size, label=label)
			elif marker == None: 	# default is a closed circle 
				ax.scatter(xsources, ysources, color=color, s=size, label=label)
		else: 
			if marker == "o": 	# 'o' used for open circle 
				plt.scatter(xsources, ysources, facecolor="none", edgecolor=color, s=size, label=label)
			elif marker == None: 	# default is a closed circle 
				plt.scatter(xsources, ysources, color=color, s=size, label=label)

	except: return "Failed to add points."		
	

###-----------------------------------------------------------------------------------------------------

def CorrectMags(frame=None, phots=None, corrections=None, field=None, apertures=[3,20], headers=[V, B, I], instrument="ACS", filters=["F606W", "F435W", "F814W"], distance=False, savefile=None, ID_header=ID, coord_headers=[X, Y], extinction=[0,0,0,0]): 

	"""Calculating magnitudes with given aperture corrections. The input 'instrument' can be 'ACS' or 'WFC3', which defines which EEF file to read from. Filters should be read in the order [V,B,I]. If given, 'extinction' should also be in [Av,Ab,Ai,Au] order. (NOTE: note RGB) or [V,B,I,U], if U is given. Corrections should also be read in VBI order. """

	try: frame = frame.copy()
	except: pass;

	if not distance: 
		distance = float(input("Distance to galaxy (in units parsec): "))

	# If U is given, add U corrections to all commands
	if len(filters) == 4: U_true = True
	else: U_true = False
	
	curr_dir = pwd()
	# Loading in the EEFs file
	cd(file_dir)
	EEFs = LoadSources(instrument + "_EEFs.txt", verbose=False)
	cd(curr_dir)

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

		Mags.to_csv(savefile)

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

def CorrectMag(df=False, phots=None, correction=None, field=None, apertures=[3,20], instrument="ACS", filt="F606W", distance=False, savefile=None, ID_header=ID, coord_headers=["X", "Y"], extinction=0): 

	"""Calculating magnitude with given aperture correction, like CorrectMags, but specifically for a single input filter (so that it doesn't require all filters to be given if only one measurement is needed). The input 'instrument' can be 'ACS' or 'WFC3', which defines which EEF file to read from. Filters should be read in the order [V,B,I]. If given, 'extinction' should also be in [Av,Ab,Ai,Au] order. (NOTE: note RGB) or [V,B,I,U], if U is given. Corrections should also be read in VBI order. """

	if df: frame = df.copy()

	curr_dir = pwd()
	# Loading in the EEFs file
	cd(file_dir)
	EEFs = LoadSources(instrument + "_EEFs.txt", verbose=False)
	cd(curr_dir)

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

	if not distance: distance = input("Input distance to galaxy (in parsecs): ")
	dmod = 5.*log10(distance/10.)

	# Converting to absolute magnitudes? 
	Mag = phots - dmod - extinction

	# Calculating the corrections
	try: 
		corr = np.array(Mag) - correction + 2.5*log10(float(EEF))
	except: 
		corr = np.array(Mag) - correction + 2.5*log10(float(EEF))

	#for i in range(len(corr)): 
	#	if corr[i] > 100 or corr[i] < -100: corr[i] = np.nan

	if savefile: 
		try: Mags = BuildFrame(headers=[ID_header, coord_headers[0], coord_headers[1], header], \
		         values = [frame[ID_header], frame[coord_headers[0]], frame[coord_headers[1]], corr])
		except: 
			Mags = BuildFrame(headers=[ID_header, coord_headers[0], coord_headers[1], header], \
		         values = [frame[ID_header], frame[coord_headers[0]], frame[coord_headers[1]], corr])

		Mags.to_csv(savefile)

	try: 
		frame[header] = corr.tolist()[0]
		return frame
	except: 
		return corr


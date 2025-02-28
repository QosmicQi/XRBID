###################################################################################
##########	For calculating and plotting hardness ratios		########### 
###################################################################################

import re
import random
import numpy as np
from numpy import pi, mean, median
import matplotlib
import matplotlib.pyplot as plt
from astropy.io.votable import parse
from numbers import Number
import pandas as pd
pd.options.mode.chained_assignment = None
from DataFrameMod import Find, FindUnique, RawFind
from Headers import ID, CSCID, T_counts, T_countslo, T_countshi, \
	U_counts, U_countslo, U_countshi, S_counts, S_countslo, \
	S_countshi, M_counts, M_countslo, M_countshi, \
	H_counts, H_countslo, H_countshi, HS, HSlo, HShi, HM, HMerr, \
	HMlo, HMhi, MS, MSerr, MSlo, MShi, \
	Flux_bb, Flux_bblo, Flux_bbhi, Flux_pow, Flux_powlo, Flux_powhi, \
	MSRat, HSRat, HMRat

imext = [0., 13500.]

# Setting the style of my plots
fontparams = {'font.family':'stix'}
labelparams = {'family':'stix', 'size':30}

###-----------------------------------------------------------------------------------------------------

def CalcHR(df, higher, lower, meanval=False, returnlabel=False): 

	"""For calculating hardness ratios of sources in a DataFrame given two hardness bands.""" 
	
	# This code assumes the user is competent and has called on the bands in the correct order. 
	# If not, the following calculations may yield odd results. 

	unique = FindUnique(df)["ID"].tolist()

	hrs = [CalcHRSub(Find(df, "ID == " + i), higher, lower, meanval) for i in unique]

	# Simplifying the label
	higher = re.split(" ", higher)[0]
	lower = re.split(" ", lower)[0]

	if returnlabel: return hrs, str("(" + higher + " - " + lower + ") / (" + higher + " + " + lower + ")")
	else: return hrs
	
###-----------------------------------------------------------------------------------------------------

def PlotHR(x0, y0, xlab=None, ylab=None, label=None, color="black", marker=None, size=10, save=None, figsize=(7,5), yscale="linear", xscale="linear", width=1, fontsize=20):

	"""Plots hardness ratio as a function of secondary input. (Can technically be used for any kind of scatterplot)"""

	x0 = np.array(x0)
	y0 = np.array(y0)

	# Removing bad values from the list, if given
	mask0 = [isinstance(i, float) for i in x0]
	mask1 = [not np.isnan(i) for i in x0[mask0]]

	y0 = y0[mask0][mask1]
	x0 = x0[mask0][mask1]

	if xlab == MSRat: xlab = "(Medium - Soft)/(Medium + Soft)"
	if xlab == HSRat: xlab = "(H - S)/(H + S)"
	if xlab == HMRat: xlab = "(Hard - Medium)/(Hard + Medium)"

	if ylab == MSRat: ylab = "(Medium - Soft)/(Medium + Soft)"
	if ylab == HSRat: ylab = "(H - S)/(H + S)"
	if ylab == HMRat: ylab = "(Hard - Medium)/(Hard + Medium)"


	plt.figure(figsize=figsize)

	if marker == "o": 	# 'o' used for open circle
		plt.scatter(x0, y0, label=label, facecolor="none", edgecolor=color, s=size**2, linewidths=width)
	elif marker == None: 	# default is a closed circle
		plt.scatter(x0, y0, label=label, color=color, s=size**2, linewidths=width)
	else: 
		plt.scatter(x0, y0, marker=marker, label=label, color=color, s=size**2, linewidths=width)

	plt.ylabel(ylab, size=20)
	plt.xlabel(xlab, size=20)
	plt.yscale(yscale)
	plt.xscale(xscale)

	plt.rcParams.update(fontparams)


	plt.tick_params(which="major", length=10, width=1.5, labelsize=fontsize, direction="inout")
	if label: plt.legend(prop={"size":18})
	if save: plt.savefig(save, dpi=500)

	#plt.show()
	

###-----------------------------------------------------------------------------------------------------

def AddHR(x0, y0, xlab=None, ylab=None, label=None, color="black", marker=None, size=10, alpha=1.0, width=1.0):

	"""Adds a plot to a pre-existing HR graph. Allows for flexible input, but X and Y should be alike between each plot."""
	
	if marker == "o": 	# 'o' used for open circle 
		plt.scatter(x0, y0, label=label, facecolor="none", edgecolor=color, s=size**2, alpha=alpha, linewidths=width)
	elif marker == None: 	# default is a closed circle 
		plt.scatter(x0, y0, label=label, color=color, s=size**2, alpha=alpha, linewidths=width)
	else: plt.scatter(x0, y0, marker=marker, label=label, color=color, s=size**2, alpha=alpha, linewidths=width)

	if xlab: plt.xlabel(xlab, size=20)
	if ylab: plt.ylabel(ylab, size=20)

###-----------------------------------------------------------------------------------------------------

def CalcHRSub(source, higher, lower, meanval=False): 
	
	"""A subroutine designed to make CalcHR more efficient. This is behind-the-scenes and should not be read imported anywhere.""" 

	temphi = Find(source, [higher + " >= 1"])[higher].tolist() # only uses sources with valid values
	templo = Find(source, [lower + " >= 1"])[lower].tolist()
	if meanval: 
		temphi = mean(temphi)
		templo = mean(templo)
	else:	# by default, use median
		temphi = median(temphi)
		templo = median(templo)
	rat = (temphi - templo)/(temphi + templo)
	if rat == 1.0: return np.nan
	else: return rat
###-----------------------------------------------------------------------------------------------------

#def 

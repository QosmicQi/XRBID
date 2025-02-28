###################################################################################
##########		For reading DataFrames for CSC sources		########### 
##########		Last update: Nov. 26, 2024 			########### 
##########		Update desc: Standardized parameters,		########### 
##########		  added descriptions, and removed unused	########### 
##########		  or obsolete code, added MakePostage()		########### 
###################################################################################

import re
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from astropy.io.votable import parse, parse_single_table
import astropy.io.votable
import glob
import pandas as pd
pd.options.mode.chained_assignment = None
import warnings
warnings.filterwarnings("ignore")
imext = [0., 13500.]
import os
cd = os.chdir
pwd = os.getcwd

from XRBID.Headers import heads, tabheads, tabheads_split, Greenx, Greeny, Redx, Redy, Bluex, Bluey, Radius, DaoNo, headers_dict, ID, Class, Class2, Conf, Notes, X, Y
from XRBID.DataFrameMod import BuildFrame, Find, Convert_to_Number

###-----------------------------------------------------------------------------------------------------

def LoadSources(infile=None, verbose=True):

	"""
	Creates and returns a DataFrame using the specified text file.
	
	PARAMETERS
	----------
	infile		[str] 		: Input file name
	verbose		[bool] (True)	: If True, prints updates on the process
	
	RETURNS
	---------
	Returns a DataFrame created from the input file
	
	"""

	if infile == None: infile = raw_input("File containing data for sources: ") 
	if len(infile.split(".")) < 2: infile = infile + ".txt";

	if verbose: print("Reading in sources from " + infile + "...") 

	try: 
		# ID will sometimes be given as a number but should be read in as string
		try: return Convert_to_Number(pd.read_csv(infile, sep=",", converters={"ID": str}, dtype=None).drop("Unnamed: 0", axis=1))
		# If ID not in the file, just return this
		except: 
			try: return Convert_to_Number(pd.read_csv(infile, sep=",", dtype=None).drop("Unnamed: 0", axis=1))
			except: print("ERROR: Reading in datafile failed.") 
	except: 
		try: return Convert_to_Number(pd.read_csv(infile, sep=",", dtype=None))
		except: 
			try: return Convert_to_Number(pd.read_csv(infile, dtype=None))
			except: 
				print("ERROR: Reading in datafile failed.") 

###-----------------------------------------------------------------------------------------------------

def NewSources(infile=None, headers=None, rename=False, outfile=False):

	"""
	Creates a new DataFrame from a specified VOTable. If headers are known, can be read in to start. 
	If rename = True, list all headers and offer user to rename in order.
	
	PARAMETERS
	----------
	infile		[str] 		: Name of the VOTable file to read from
	headers		[list]		: Optional; list of header names
	rename		[bool] (False)	: If True, allows user to rename the headers
	outfile		[str] (False)	: Name of save file. 
	
	RETURNS
	---------
	sources		[pd.DataFrame]	: DataFrame created from the input file

	""" 

	# user must give VOTable file to read from
	if infile == None: infile = input("VOTable file: ") 
	if len(infile.split(".")) < 2: infile = infile + ".vot"

	print("Reading in table from " + infile)

	votable = parse(infile)   # reads in table 
	table = votable.get_first_table().to_table(use_names_over_ids=True)

	values = []		# Keeps track of the values in the votable

	# Keeping track of the headers from votable
	temp_headers = [f.name for f in parse_single_table(infile).fields]

	# Populating values and temp_headers
	for i in temp_headers: 
		values.append(table[i])


	# If no headers are speficied, read in all in the votable file
	if not headers and rename:
		print("Enter new header names, or hit 'Enter' to keep name.\n")
		frame_heads = []
		for i in temp_headers:
			temp_head = input(i+": ")
			if not temp_head: 
				temp_head = i
			frame_heads.append(temp_head)
	elif headers: frame_heads = headers
	else: frame_heads = temp_headers

	sources = BuildFrame(headers=frame_heads, values=values)

	# For some reason, reading directly from a CSC .vot does not allow
	# Find() to be used on the DataFrame unless read back in by LoadSources.
	# Saving to a file, then reading back in with LoadSources. 
	if outfile: 
		sources.to_csv(outfile)
		sources = LoadSources(outfile)
	else: 
		sources.to_csv("temp.txt")
		sources = LoadSources("temp.txt", verbose=False)
		os.remove("temp.txt")

	print("DONE")
	return sources

###-----------------------------------------------------------------------------------------------------

def GetCoords(infile=None, IDs=None, savecoords=None, checkcoords=False, verbose=True):

	"""
	Gets the image (X/Y) coordinates of the sources from a given region, text file, or DataFrame file. 
	The IDs argument can be used to only return certain coordinates, such as in the case of finding DaoFind 
	source coordinates from a region file of specific DaoFind source IDs. 
	NOTE: to enable this to work in non-Daofind cases, will need to edit code some more.
	
	PARAMETERS
	----------
	infile		[str] 		: Name of the VOTable file to read from
	headers		[list]		: Optional; list of header names
	rename		[bool] (False)	: If True, allows user to rename the headers
	
	RETURNS
	---------
	sources		[pd.DataFrame]	: DataFrame created from the input file

	"""

	if verbose: print("Retrieving coordinates from " + infile)

	if ".frame" in infile: # DataFrame files can be read in using LoadSources
		temp = LoadSources(infile=infile, verbose=False)
		try: 	# Let image coordinates be the default
			x_coords = temp[X].values
			y_coords = temp[Y].values
		except: # If image coordinates not available, use RA and Dec.
			x_coords = temp[RA].values
			y_coords = temp[Dec].values
	else: 
		try: 	# If the file is a .txt file (or some other), this should work fine
			coords = np.genfromtxt(infile)
			if checkcoords: 
				userin = "n"
				print("Input check on " + infile + ": ")
				i = -1
			else: 
				userin = "y"
				i = 0
			
			while userin != "y": 
				i += 1
				# Trying to find the coordinates in the file. 
				# This should cover all kind of file formats
				try: 
					userin = raw_input("Is (" + str(coords[0][i]) + ", " + str(coords[0][i+1]) + ") a valid coordinate? ").lower()[0]
				except: print("No coordinates found. Check file."); break;
				
			x_coords = coords.T[i]
			y_coords = coords.T[i+1]
		except:	# If the file is a region file (.reg), the following should apply 
			with open(infile) as f: 
			    lines = f.readlines()[3:]
			    
			x_coords = []
			y_coords = []

			for i in range(len(lines)):
				line = lines[i]
				try: 
					l0 = float(line.strip().split()[0].split('(')[-1].split(',')[0])
					l1 = float(line.strip().split()[1].split(",")[0])
				except: 
					l = line.strip().split()[0].split('(')[-1].split(',')
					l0 = float(l[0])
					l1 = float(l[1])
				x_coords.append(l0)
				y_coords.append(l1)
	
	if IDs != None: 
		# Create a mask of all False. Then replace the corresponding ID with True to mask out all coordinates that do not correspond with the sources in IDs.
		mask = [False]*len(x_coords)
		for i in IDs: 
			j = int(i) - 1
			mask[j] = True
		# Applying the mask
		x_coords = np.array(x_coords)[mask]
		y_coords = np.array(y_coords)[mask]
		x_coords = x_coords.tolist()
		y_coords = y_coords.tolist()

	if savecoords: 
		print("Saving " + savecoords)
		with open(savecoords, "w") as f: 
			np.savetxt(f, np.column_stack([x_coords, y_coords]))
	return x_coords, y_coords

###-----------------------------------------------------------------------------------------------------

def GetIDs(infile=None, verbose=True):

	"""
	Gets the IDs (or printed text) of the sources from a given region or text file.
	
	PARAMETERS
	----------
	infile	[str]		: Name of input file
	verbose	[bool] (True)	: If True, prints file info

	RETURNS
	---------
	ids	[list]		: List of source IDs
	
	"""

	if verbose: print("Retrieving IDs from " + infile)

	# If the file is a region file (.reg), the following should apply 
	with open(infile) as f: 
	    lines = f.readlines()[3:]
	    
	ids = []

	for i in range(len(lines)):
		line = lines[i]
		try: 
			ids.append(re.split("}", re.split("text={", line)[-1])[0])
		except: 
			ids.append("None")

	return ids

###-----------------------------------------------------------------------------------------------------

def SourceID(infile=None): 

	"""
	Reads in initial identifications from the *_IDs.txt file and returns a new DataFrame. 
	Once this is read in once, the file may be saved as a .txt file that can be read by LoadSources().
	This function has a very specific use and is not useful for most people. Please ignore!
	"""

	if infile == None: infile = raw_input("ID file name? : ")
	if len(infile.split(".")) <= 1: infile = infile + ".txt" 
	
	try: 
		with open(infile) as f: 
			lines = f.readlines()[1:]
	except: print("File not found.") 

	sourceid = pd.DataFrame(columns=[ID, DaoNo, Class, Class2, Conf, Notes])
	sourceid[ID] = [" "]*len(lines) # setting the length of the DataFrame

	try: 
		try: 
			for i in range(len(lines)): 
				j = lines[i].split()
				sourceid[ID][i] = j[0]
				sourceid[DaoNo][i] = j[1]
				sourceid[Class][i] = j[2]
				sourceid[Class2][i] = j[3]
				sourceid[Conf][i] = float(j[4])
				try: sourceid[Notes][i] = " ".join(j[5:])
				except: sourceid[Notes][i] = "None"
		except: print("Error on " + j[0] + ". DataFrame not created.")

		for i in range(len(sourceid)): 
			# Resetting classifications to something more readable
			if sourceid[Class][i] == "H": sourceid[Class][i] = "HMXB"
			elif sourceid[Class][i] == "I": sourceid[Class][i] = "IMXB"
			elif sourceid[Class][i] == "L": sourceid[Class][i] = "LMXB"
			elif sourceid[Class][i] == "A": sourceid[Class][i] = "AGN"
			elif sourceid[Class][i] == "S": sourceid[Class][i] = "SNR"
			elif sourceid[Class][i] == "Q": sourceid[Class][i] = "Quasar"
			elif sourceid[Class][i] == "C": sourceid[Class][i] = "Cluster"
			elif sourceid[Class][i] == "?": sourceid[Class][i] = "Unknown" 

			if sourceid[Class2][i] == "H": sourceid[Class2][i] = "HMXB"
			elif sourceid[Class2][i] == "I": sourceid[Class2][i] = "IMXB"
			elif sourceid[Class2][i] == "L": sourceid[Class2][i] = "LMXB"
			elif sourceid[Class2][i] == "A": sourceid[Class2][i] = "AGN"
			elif sourceid[Class2][i] == "S": sourceid[Class2][i] = "SNR"
			elif sourceid[Class2][i] == "Q": sourceid[Class2][i] = "Quasar"
			elif sourceid[Class2][i] == "C": sourceid[Class2][i] = "Cluster"
			elif sourceid[Class2][i] == "X": sourceid[Class2][i] = "N/A"
			elif sourceid[Class2][i] == "?": sourceid[Class2][i] = "Unknown"

		return sourceid

	except: print("Error creating DataFrame.")

###-----------------------------------------------------------------------------------------------------


def SourceList(savefile, df=None, columns=['ID']):

	""" 
	Creates a text file containing some input data stacked into columns. 
	By default, if a DataFrame is read in and no columns are specified, 
	the code searches for a header in the DataFrame called 'ID' and prints
	those values to a file. Otherwise, SourceList() can be used to print
	the values of any header within DataFrame given in a list called 
	columns, or columns can be used as a list of values to print to an
	output file (if not DataFrame is given).
	(Modified for simplification Nov 26, 2024)
	
	PARAMETERS
	----------
	savefile	[str]		: Name of output file
	df		[pd.DataFrame] 	: DataFrame containing information to save. Optional only if columns is given.
	columns		[list] (['ID'])	: Optional; contains either the names of the headers to pull from the input
					  DataFrame, or a list of values to save to the output file. 
	
	RETURNS
	---------
	Saves a file under the name savefile

	""" 

	if df: 
		with open(savefile, 'w') as f:
			tempstack = []
			for i in columns: 
				tempstack.append(df[i].values)
			np.savetxt(f, np.column_stack(tempstack), fmt="%s")
	else: np.savetxt(f, np.column_stack(columns), fmt="%s")

	print(savefile + " saved!")

###-----------------------------------------------------------------------------------------------------

def DaoClean(daosources=None, sources=None, wiggle=0): 

	"""
	Cleaning the DaoFind sources to exclude any candidate that falls outside of the radius of the 
	X-ray sources, taking into account given wiggle room (in pixels). Sources should be read in as
	a dataframe with image coordinates (X, Y) and 2-sig radius saved under header Radius.
	
	PARAMETERS
	----------
	daosources	[df.DataFrame]	: DataFrame containing the coordinates of sources identified DaoFind,
					  where the coordinates should be saved under the headers 'X' and 'Y'
	sources		[df.DataFrame]	: DataFrame containing the X-ray sources, with 2sigma radii saved under
					  the header 'Radius' and coordinates under 'X' and 'Y'
	wiggle		[float] (0)	: Additional pixels to add to the search radius, giving a more lenient search
	
	RETURNS
	---------
	GoodSources	[df.DataFrame]	: DataFrame equivalent to 'daosources' with only sources that fall within
					  the radii of the sources in the DataFrame 'sources.'

	"""


	try: daosources = daosources.copy()
	except: pass;

	# Retrieving headers from the DataFrames
	daoheads = daosources.columns.values.tolist()
	headlist = daosources.columns.values.tolist()
	headlist.append("ID")

	### Cleaning daofind sources to only include those around our sample sources
	daocleaned = np.empty((0,len(headlist)))

	print("Cleaning DAOFind sources. This will take a few minutes. Please wait.. ")
	for i in range(len(sources)): 
		# Properties of each Lehmer source
		temprad = sources[Radius][i]
		xtemp = sources[X][i]
		ytemp = sources[Y][i]
		tempid = sources[ID][i]

		if xtemp >=0 and ytemp >= 0:
			# Search area around each source
			#print(temprad, xtemp, ytemp)
			tempxmax = xtemp+temprad+wiggle
			tempxmin = xtemp-temprad-wiggle
			tempymax = ytemp+temprad+wiggle
			tempymin = ytemp-temprad-wiggle

		#print(tempxmax, tempxmin, tempymax, tempymin)
		# Finding the coordinates of daofind sources within the square area of the Lehmer source
		try: 
			tempdao = Find(daosources, ["green x >= " + str(tempxmin), "green x <= " + str(tempxmax), \
				                 "green y >= " + str(tempymin), "green y <= " + str(tempymax)])

			for j in range(len(tempdao)): 
				if sqrt((tempdao[Greenx][j] - xtemp)**2 + (tempdao[Greeny][j] - ytemp)**2) <= temprad+wiggle:
					tempstack = [tempdao[k][j] for k in daoheads]
					tempstack.append(tempid)
					daocleaned = np.vstack((daocleaned, tempstack))
		except: 
			tempdao = Find(daosources, ["x >= " + str(tempxmin), "x <= " + str(tempxmax), \
				                 "y >= " + str(tempymin), "y <= " + str(tempymax)])

			for j in range(len(tempdao)): 
				if sqrt((tempdao[X][j] - xtemp)**2 + (tempdao[Y][j] - ytemp)**2) <= temprad+wiggle:
					tempstack = [tempdao[k][j] for k in daoheads]
					tempstack.append(tempid)
					daocleaned = np.vstack((daocleaned, tempstack))

	print("DONE WITH CLEANING. CREATING DATAFRAME...")

	daotemp = [daocleaned.T[k] for k in range(len(daocleaned.T))]
	GoodSources = BuildFrame(headers=headlist, values=daotemp)

	return GoodSources

###-----------------------------------------------------------------------------------------------------

def GetDaoPhots(df=None, ann="10-3", gal="M83", getcoords=True): 

	"""
	Retreiving the photometry for each Dao No source in a given DataFrame. 
	Returns V, B, I, V-I, B-V, and B-I photometry (already adjusted) from 
	phots_f#_ann##-#.frame files. 	
	
	Example use: 
	A DataFrame is created containing the optical photometry of all DaoFind sources in HST field 1,
	with a background annulus of radius 3 to 10 pixels, and saved as phots_f1_ann10-3.frame.
	DaoClean() returns a second DataFrame with the IDs (Dao No) of candidate optical counterparts. 
	Running GetDaoPhots on the second DataFrame 

	PARAMETERS
	----------
	df	[df.DataFrame]	: DataFrame containing the coordinates of sources identified DaoFind,
					  where the coordinates should be saved under the headers 'X' and 'Y'
	sources		[df.DataFrame]	: DataFrame containing the X-ray sources, with 2sigma radii saved under
					  the header 'Radius' and coordinates under 'X' and 'Y'
	wiggle		[float] (0)	: Additional pixels to add to the search radius, giving a more lenient search
	
	RETURNS
	---------
	GoodSources	[df.DataFrame]	: DataFrame equivalent to 'daosources' with only sources that fall within
					  the radii of the sources in the DataFrame 'sources.'

	
	"""

	sources = df.copy()

	print("Loading Daofind source magnitudes. May take a few minutes.")
	
	# Loading the sources for each field

	Dao1 = LoadSources("daomags_f1_ann"+ann+".frame")
	Dao2 = LoadSources("daomags_f2_ann"+ann+".frame")
	Dao3 = LoadSources("daomags_f3_ann"+ann+".frame")
	Dao4 = LoadSources("daomags_f4_ann"+ann+".frame")
	Dao5 = LoadSources("daomags_f5_ann"+ann+".frame")
	Dao6 = LoadSources("daomags_f6_ann"+ann+".frame")
	Dao7 = LoadSources("daomags_f7_ann"+ann+".frame")

	# Creating a searchable list of daosource DataFrames
	Daos = [Dao1, Dao2, Dao3, Dao4, Dao5, Dao6, Dao7]

	# List for gathering necessary values
	v = []
	b = []
	i = []
	vi = []
	bv = []
	bi = []
 
	if getcoords: 
		Coords1 = GetCoords(gal+"_f1_f555w_daosources.reg") 
		Coords2 = GetCoords(gal+"_f2_f547m_daosources.reg") 
		Coords3 = GetCoords(gal+"_f3_f547m_daosources.reg") 
		Coords4 = GetCoords(gal+"_f4_f547m_daosources.reg") 
		Coords5 = GetCoords(gal+"_f5_f547m_daosources.reg") 
		Coords6 = GetCoords(gal+"_f6_f547m_daosources.reg") 
		Coords7 = GetCoords(gal+"_f7_f547m_daosources.reg")
		Coords = [Coords1, Coords2, Coords3, Coords4, Coords5, Coords6, Coords7]
		xcoords = []
		ycoords = []

	print("\nGetting source magnitudes...\n")

	# Going through each of the sources in the given DataFrame
	for j in range(len(sources)):
		
		# Gathering the field and Dao No. of current source
		f = str(sources["Field"][j])
		dao = sources["Dao No"][j]


		# Choosing the correct Daosources using the field number
		# Daos is a list of the Daomags
		DaoTemp = Daos[int(f) - 1]
		
		# Choosing the sourcce from Dao No. (daomags does not have Dao No.)
		# DaoTemp = Find(Dao, "Dao No = " + dao)
		
		# Gathering and calculating magnitudes
		try: 
			v.append(DaoTemp["V"][dao])
			b.append(DaoTemp["B"][dao])
			i.append(DaoTemp["I"][dao])
			vi.append(v[-1] - i[-1])
			bv.append(b[-1] - v[-1])
			bi.append(b[-1] - i[-1])
			# For some reason, this is wrong in the file, and I don't understand why
			#vi.append(DaoTemp["V-I"][dao])
			#bv.append(DaoTemp["B-V"][dao])
			#bi.append(DaoTemp["B-I"][dao])
		except: # If no sources is given (in case of LMXB), set as NaN
			v.append(np.nan)
			b.append(np.nan)
			i.append(np.nan)
			vi.append(np.nan)
			bv.append(np.nan)
			bi.append(np.nan)

		# Gather coords if requested
		if getcoords: 
			try:
				xcoords.append(Coords[int(f)-1][0][int(dao)])
				ycoords.append(Coords[int(f)-1][1][int(dao)])
			except: 
				xcoords.append(np.nan)
				ycoords.append(np.nan)

	print("Done!")
	if getcoords: 
		return v, b, i, vi, bv, bi, xcoords, ycoords
	else: return v, b, i, vi, bv, bi

###-----------------------------------------------------------------------------------------------------

def ListFiles(search='*'): 

	""" 
	For getting a list of files from the current directory.

	PARAMETERS
	----------
	search	[str] ('*')	:  Search criteria to narrow down the list of files
	
	RETURNS
	---------
	temp 	[list]		: List of files matching the search criteria

	"""

	temp = []
	for file in glob.glob(search): 
		temp.append(file)
	temp.sort()
	return temp

###-----------------------------------------------------------------------------------------------------

def Crossref(df=None, regions=None, coords=None, outfile="crossref_results.txt", search_radius=30, catalogs=None, verbose=True, getcoords=False, coordsys="img", coordheads=None): 

	"""
	From input DataFrame and/or region files (in image coordinate format), finds overlaps within a given 
	search radius of the DataFrame sources and prints all ID names to a file as a DataFrame. 
	If the coordinates are given as [RA, Dec] instead of [X,Y], must change coordsys from "img" to "fk5" 
	and convert search_radius from pixels to degrees. Can feed in the name of the catalogs used to output 
	as DataFrame headers. Otherwise, the region name will be used.

	NOTE: There is an error in this where if the first region file doesn't have a counterpart in the first 
	entry of the overlap file, the first entry may be split into multiple entries. Check file.
	
	"""

	sources = df.copy()

	xlist = []
	ylist = []
	idlist = []

	if not isinstance(search_radius, list): search_radius = [search_radius]*len(sources)

	masterlist = [] # list of all matches sources

	if regions:
		if not isinstance(regions, list): regions = [regions]
		for i in regions: 
			idlist.append(GetIDs(i, verbose=False))
			xtemp, ytemp = GetCoords(i, verbose=False)
			xlist.append(xtemp)
			ylist.append(ytemp)
	elif coords: 
		# if given coords, they should be read in as a list of [xcoords, ycoords]. 
		xlist = coords[0]
		ylist = coords[1]
		if not isinstance(xlist, list): 
			xlist = [xlist]
			ylist = [ylist]

	blockend = 0
	if verbose: print("Finding cross-references between sources. This will take a few minutes. Please wait.. ")
	for i in range(len(sources)): # for each source in the DataFrame
		# Properties of each Lehmer source
		if coordsys == "img" and coordheads == None: 
			try: 
				xtemp = sources[X][i]
				ytemp = sources[Y][i]
			except: 
				xtemp = sources["x (mosaic)"][i]
				ytemp = sources["y (mosaic)"][i]
		elif coordsys == "fk5" and coordheads == None: 
			xtemp = sources["RA"][i]
			ytemp = sources["Dec"][i]
		elif coordheads: 
			xtemp = sources[coordheads[0]]
			ytemp = sources[coordheads[1]]

		tempid = sources[ID][i]
		tempn = 0  
		# tempn keeps track of the number of overlap sources identified in the current list for the current base source (used as index)

		#if xtemp >=0 and ytemp >= 0: # Removed because I want to print all crossrefs, not just 'valid' sources within the image.
		# Search area around each source
		tempxmax = xtemp+search_radius[i]
		tempxmin = xtemp-search_radius[i]
		tempymax = ytemp+search_radius[i]
		tempymin = ytemp-search_radius[i]

		#print(tempxmax, tempxmin, tempymax, tempymin)

		# Adding each new source to the list
		# If no counterparts are found, the source will appear with "None" for counterparts.
		tempids = [None]*(len(idlist) + 1)
		tempids[0] = tempid 
		masterlist.append(tempids)

		# Searching each list of sources from each region file to identify overlaps
		for j in range(len(idlist)): # Number of lists (region files) to search through 
			for k in range(len(xlist[j])): # Number of sources to search through for the current list/region file
				# When overlap is found, see if masterlist has room to add it. 
				# If not, add a new row to make room.
				if tempxmax > xlist[j][k] > tempxmin and tempymax > ylist[j][k] > tempymin and \
				sqrt((xlist[j][k]-xtemp)**2 + (ylist[j][k]-ytemp)**2) <= search_radius[i]: 
					try: 
						# With blockend showing how many total items were found prior to the search on this source, 
						# and tempn showing how many counterparts were identified or the current source, 
						# blockend+tempn should identify the index of the current source
						# The following will cycle through all indices from blockend to blockend+tempn 
						# to see where the last open space is
						for n in range(tempn+1):
							if masterlist[blockend+n][j+1] == None: 
								masterlist[blockend+n][j+1] = idlist[j][k]
								break; # After the last open space, break the chain.
							else: pass;
					except: 
						# Exception will be raised once we reach the end of the current list without finding a free space. 
						# Add a new line, if that's the case.
						tempids = [None]*(len(idlist) + 1) # keeps track of the ids associated with the identified source
						tempids[0] = tempid
						tempids[j+1] = idlist[j][k]
						masterlist.append(tempids)

					tempn = tempn + 1 # adds a count to the identified sources for this file.

		blockend = len(masterlist) # the index of the end of the previous "block" of sources already detected.

	if verbose: print("DONE WITH CLEANING. CREATING DATAFRAME...")

	# If catalogs not given, use the name of the region files as the headers of the DataFrame
	if not catalogs:
		catalogs = []
		for i in regions: catalogs.append(i.split(".reg")[0]) 

	# Adding catalogs to the headers to be read into the DataFrame
	headlist = [ID]
	for i in catalogs: headlist.append(i)

	vallist = []
	# Converting the masterlist into an array to be read in as DataFrame values
	temp_array = np.array(masterlist).T
	for i in range(len(temp_array)): vallist.append(temp_array[i].tolist())

	GoodSources = BuildFrame(headers=headlist, values=vallist)
	GoodSources.to_csv(outfile)

	return GoodSources

###-----------------------------------------------------------------------------------------------------

def GalComponents(sources, rad=[0], locs=["Disk", "Outskirt"], theta=0, center=None, savereg=False, regname="sources"): 

	"""
	Breaks up the sources within a given galaxy into the the regional components of the galaxy as given by the 'locs' argument. 
	The radii should be given in pixel units. It may be read in as a list of radii, in order of innermost regions outward. 
	The argument 'locs' should give the corresponding location names and should be one element larger than the size of the radius 
	list (preferrably with outskirt being the last element, corresponding to all sources outside of the defined regions of interest). 
	If a given element in the rad list is a list of multiple radii, assume an ellipse at an angle theta and run the ellipse 
	function InEllipse. GalComponents returns the same DataFrame with an added header, Location, which details which component 
	the source appears within.


	"""


	from WriteScript import WriteReg

	sources = sources.copy()
	locations = ["Outskirt"]*len(sources)
	
	# Finding the center of the galaxy using either given coordinates or the nucleus of the galaxy.
	if center: 
		xcenter = center[0]
		ycenter = center[1]
	else: 
		temp = Find(sources, "Class = Nucleus")
		try: 
			xcenter = temp[X][0]
			ycenter = temp[Y][0]
		except: 
			xcenter = temp["x (mosaic)"][0]
			ycenter = temp["y (mosaic)"][0]

	## Start with the outermost region and work inwords
	for i in range(len(sources)): 
		try: 
			xtemp = sources[X][i]
			ytemp = sources[Y][i]
		except: 
			xtemp = sources["x (mosaic)"][i]
			ytemp = sources["y (mosaic)"][i]

		#starting with the outermost region and works way inword
		for j in range(1, len(rad)+1): 
			r = rad[-j] # first element is the last in the radius list
			if not isinstance(r, list): # if only one radius given, assume circle
				if sqrt((xcenter - xtemp)**2 + (ycenter - ytemp)**2) <= r:
					# ignoring the last element of the locs list, save the location
					# of the next outermost region
					locations[i] = locs[-j-1]
				else: pass;
			else: # if multiple radii are given, assume ellipse 
				if InEllipse(xtemp, ytemp, [xcenter, ycenter], r[0], r[1], theta): 
					locations[i] = locs[-j-1]
				else: pass;

	sources["Location"] = locations

	if savereg: 
		for i in locs: 

			print(i)
			WriteReg(Find(sources, "Location = " + i), outfile=regname+"_"+i.lower()+".reg", color="red", width=2, radius=50)		

	print("Locations of sources successfully identified")
	return sources

###-----------------------------------------------------------------------------------------------------

def CheckBounds(df, imext=imext, remove=False, search=["x", "y"], resetbounds=False, ID=ID):

	""" 
	Checking whether sources in the given DataFrame is inside the bounds in the image. 
	Returns a DataFrame with the boundary conditions (whether the source is in or out) marked. 

	""" 

	# For now, assume Bounds is given. If it isn't, set remove to False and only return sources within the bounds.

	frame = df.copy()

	try: 
		if resetbounds: frame[Bounds] = "Out" # resets bounds
	except: remove = True

	# converting image extent (imext) to an array
	# if the length of the array is 1, use imext[0] as min and imext[1] as max for both x and y
	# otherwise, the length should 2, for x-extent and y-extent
	imext = np.array(imext)
	if len(imext) == 1: 
		xmin = ymin = imext[0]
		xmax = ymax = imext[1]
	else: 
		xmin = imext[0][0]
		xmax = imext[0][1]
		ymin = imext[1][0]
		ymax = imext[1][1]
	
	# finding all of the sources within bounds
	temp = Find(frame, [search[0] + " >= " + str(xmin),  search[0] + " <= " + str(xmax), search[1] + " >= " + \
		    str(ymin), search[1] + " <= " + str(ymax)])

	if remove: 
		frame = temp
		try: 
			for i in range(len(frame)): 
				frame[Bounds][i] = "In"
		except: pass; 
	else:	# if not removing all out of bounds, then find all 
		try: inID = temp[ID].values.tolist()
		except: ID = raw_input("ID not found. Enter ID Header: ") 

		inID = temp[ID].values.tolist()  # list of in-bound sources
		frame = frame.set_index(frame[ID].values) # Allows frame to be searchable by ID name

		# setting all in-bound sources to "In"
		for i in inID: 
			frame[Bounds][i] = "In" 

	# resetting the indices of the frame to source number
	frame = frame.set_index(np.arange(len(frame)))
	print("DONE")
	return frame

###-----------------------------------------------------------------------------------------------------

def InEllipse(x,y,center,rad1,rad2,theta): 

	"""Checks whether a source falls within a given ellipse and returns a bool."""

	# Calculating the position of the source relative to the center of the galaxy
	# If pos =< 1, the source is within the ellipse
	theta = theta * 0.0174533
	pos = (((x-center[0])*np.cos(theta)+(y-center[1])*np.sin(theta))**2)/rad1**2 + \
	      (((x-center[0])*np.sin(theta)-(y-center[1])*np.cos(theta))**2)/rad2**2

	# Return True with the source is within the ellipse
	return pos <= 1

###-----------------------------------------------------------------------------------------------------

def FindField(df=None, coords=None, fieldframes=["M81_fields_red.frame", "M81_fields_green.frame", "M81_fields_blue.frame"]):

	""" 
	Automatically detects the best frames that each X-ray source falls within, for red, green, and blue. Requires you 
	to have already created a Dataframe file for each filter of the fields which gives the minimum and maximum of the 
	x and y coordinates of each field. Assumes the fields are a perfect square, so will need to be double-checked. 
	Also assumes all coordinates are in pixels (img). Currently, this code only finds the coordinates for three fields. 
	Will need to be modified if more are needed! In the case where the source falls within more than one field, it takes 
	the first field of the red or the one nearest the red field for green and blue. Currently returns a list of the fields 
	(will need to modify to return a DataFrame instead).

	"""

	if len(df) > 0: coords = [df["x"].values.tolist(),df["y"].values.tolist()]
	else: pass;
	
	redfield = []
	greenfield = []
	bluefield = []

	# Loading the coordinates of the fields in each filter
	FieldsRed = LoadSources(fieldframes[0])
	FieldsGreen = LoadSources(fieldframes[1])
	FieldsBlue = LoadSources(fieldframes[2])

	for i in range(len(coords[0])): 
		tempx = coords[0][i]
		tempy = coords[1][i]
		
		# Finding which field source falls within
		Temp = Find(FieldsRed, ["Min x <= " + str(tempx), "Max x >= " + str(tempx), "Min y <= " + str(tempy), "Max y >= " + str(tempy)])

		# For the red field, choose the first field
		if len(Temp) == 0: redfield.append(np.nan)
		else: redfield.append(Temp["Field"][0])
	
		# Finding the Green frame
		Temp = Find(FieldsGreen, ["Min x <= " + str(tempx), "Max x >= " + str(tempx), "Min y <= " + str(tempy), "Max y >= " + str(tempy)])
		
		# If none found, set to np.nan. 
		# If more than one is found, default to the one that matches redfield or the first, 
		# if there is no red field given

		if len(Temp) == 0: greenfield.append(np.nan)
		else: 
			Temp2  = Find(Temp, "Field = " + str(redfield[-1]))
			if len(Temp2) == 0: greenfield.append(Temp["Field"][0])
			else: greenfield.append(Temp2["Field"][0])

		# Finding the Blue frame
		Temp = Find(FieldsBlue, ["Min x <= " + str(tempx), "Max x >= " + str(tempx), "Min y <= " + str(tempy), "Max y >= " + str(tempy)])		

		# If none found, set to np.nan. 
		# If more than one is found, default to the one that matches redfield or the first, 
		# if there is no red field given

		if len(Temp) == 0: bluefield.append(np.nan)
		else: 
			Temp2  = Find(Temp, "Field = " + str(redfield[-1]))
			if len(Temp2) == 0: bluefield.append(Temp["Field"][0])
			else: bluefield.append(Temp2["Field"][0])

	print("Finished. Please manually verify the field assignments.")

	# Return the list of fields. Remember to double-check!
	return [redfield,greenfield,bluefield]

###-----------------------------------------------------------------------------------------------------

def Mosaic(df=None, sourceimgs=None, findimgs=None, rows=6, columns=5, filename="mosaic_temp", toplabels=None, bottomlabels=None, top_coords=[225,60], bottom_coords=[225,420], remove_env=True, fontsize=20):

	"""
	Plots a mosaic of the given sources. Either a DataFrame of the sources or a list of the images can be given. Alternatively, 
	can read in a search term [e.g. *.png] to find the appropriate files. 
	NOTE: I don't think the search works currently, probably because of the last step of the search process (4/25/22).

	"""

	# If sourceimgs not given, search for the files with fingimgs	
	if findimgs: 
		sourceimgs = []
		for file in glob.glob(findimgs): 
			sourceimgs.append(file)
		sourceimgs.sort()
	
		# If a DataFrame is given, only use the images of sources within that DataFrame
		tempimgs = []
		for i in sourceimgs: 
			temp = i.split(findimgs.split("*")[0])[1].split(".")[0].split("_")[0] # Finds the ID of the image
			# This assumes the name of the file puts the ID right after *
			if len(Find(df, "ID = " + str(temp))) > 0: tempimgs.append(i)
		sourceimgs = tempimgs

	# Remove the environmental snapshots from the sourceimg list.
	for i in sourceimgs: 
		if remove_env and "env" in i: sourceimgs.remove(i)

	# Calculates the number of mosaics needed to contain all images in the given format
	addone = len(sourceimgs) % (rows*columns) > 0
	totmos = round(len(sourceimgs)/(rows*columns)) + addone

	# If labels isn't a list, assume it's a header.
	if not isinstance(toplabels, list): 
		try: toplabels = df[toplabels].values.tolist()
		except: toplabels=None


	if not isinstance(bottomlabels, list): 
		try: bottomlabels = df[bottomlabels].values.tolist()
		except: bottomlabels=None

	i = 0
	for l in range(1,totmos+1):
		f, ax = plt.subplots(rows,columns, figsize=(columns*3-1,rows*3-1))
		for j in range(0,rows):
			
			for k in range(0,columns):
				try:
					ax[j,k].axis("off")
					ax[j,k].set_aspect("equal")
					ax[j,k].imshow(plt.imread(sourceimgs[i]))
					if toplabels: 
						txttop = ax[j,k].text(top_coords[0], top_coords[1], toplabels[i], color='white', ha='center', weight="extra bold", size=fontsize)
						txttop.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
					if bottomlabels: 
						txtbot = ax[j,k].text(bottom_coords[0], bottom_coords[1], bottomlabels[i], color='white', ha='center', weight="extra bold", size=fontsize)
						txtbot.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
				except: pass;
				i += 1
		plt.subplots_adjust(wspace=0, hspace=0.02)
		plt.savefig(filename+"_"+str(l)+".png", dpi=300, bbox_inches="tight")
		plt.show()
	
	print (filename+"_*.png saved!")


###########################################################################################
##########	For writing scripts, such as bash scripts and region files	        ########### 
##########	                Last Update: Feb 11, 2025               	        ########### 
##########	(Standardized parameters, added descriptions, cleaned up code)      ########### 
###########################################################################################

import datetime
import re
import numpy as np
from astropy.io.votable import parse
import pandas as pd
pd.options.mode.chained_assignment = None
from XRBID.DataFrameMod import FindUnique
import math
from XRBID.Sources import GetCoords
from XRBID.DataFrameMod import BuildFrame

X = "x" 
Y = "y"
RA = "RA" 
Dec = "Dec" 
ID = "ID"

pixtoarcs = 0.5/12.626263
arcstodeg = 0.000277778  # 1 arcs = 0.000277778 deg
pixtodeg = arcstodeg * pixtoarcs

#print("Load test")

###-----------------------------------------------------------------------------------------------------

def WriteDS9(frame=None, galaxy="M81", colorfiles=None, regions=None, scales="zscale", imgnames=None, imgsize=[305,305], outfile=None, unique_scale=False, coords=None, ids=None, idheader="ID", filetype="jpeg", basefilter=["red"], filterorder=["red", "green", "blue"], zoom=8, env_zoom=2, coordsys=None, tile=False): 

	"""Arguments: (DataFrame, Galaxy name, Galaxy color .fits files in RGB order, Region files, Desired color scales, Output image name (number will be appended to the image), Output file name).

 Generates a bash script that will write a program for taking thumbnails/snapshots in DS9 of sources in a given DataFrame. Will take in one or more region files to open over the given galaxy file, and different scales to use to adjust the RGB colors of the image. Scales should be input with the format ["scale name", redscale, greenscale, bluescale], or, if unique_scale is True, the name of the file with the unique scalings should be used. If zscale is being used, only ["zscale"] is needed (**** SETTING TO ZSCALE IS CURRENTLY BROKEN. INSTEAD, ALWAYS USE A UNIQUE SCALES***). Parameter basefilter can be used to set which filter the region file should be aligned to, if not the default green (Ex. red used as the base filter in M81). The very first basefilter should be the main filter used for imaging, even if the region files are aligned to different filters. Filters should be called in the order they're intended to be used (default [red, green, blue]. If a different order is used, the order should be declared using the 'filterorder' parameter. """
	

	if ".txt" in scales: unique_scale = True

	if tile: # If tiling is set to True, need to run this enitrely differently
		WriteDS9_tile(frame=frame, galaxy=galaxy, colorfiles=colorfiles, regions=regions, scales=scales, \
			      imgnames=imgnames, imgsize=imgsize, outfile=outfile, unique_scale=unique_scale, \
		              coords=coords, filetype=filetype, zoom=zoom, env_zoom=env_zoom, coordsys=coordsys,ids=ids)

	# If tiling is not set to true, continue...
	else: 
		good = False

		# If only coordinates are given, build DataFrame from those.
		try: 
			if not frame: 
				if coordsys == "fk5" or coordsys == "wcs": 
					headers = ["ID", "RA", "Dec"]
				else: headers = ["ID", "x", "y"]
				if not ids: ids = np.arange(coords[0]).tolist()
				frame = BuildFrame(headers=headers, values=[ids, coords[0], coords[1]])
		except: pass;
	
		frame = frame.copy()
		display(frame) 

		if colorfiles == None: # Name of the files for RGB frames needed
			while not good: 
				colorfiles = input("Names of .fits files for galaxy color map (red, green, blue; separate by commas): ")
				colorfiles = colorfiles.split(",")
				if len(colorfiles) == 3: good = True
				else: print("Not enough color files provided."); colorfiles = None; 
		
		if regions == None: # In case region files were forgotten, ask user for input. 
			regions = input("Region files? (separate by commas, or Enter if none):") 
			regions = regions.split(",")
		if not isinstance(regions, list): regions = [regions]  # region files should be lists

		# basefilter should be indicated for every region
		# If only one basefilter is given, assume all regions are aligned to that filter.
		if not isinstance(basefilter, list): 
			basefilter = [basefilter]
		if len(basefilter) == 1 and regions: basefilter = [basefilter[0]]*len(regions)
		base = filterorder[0]

		#if imgnames == None: # Names of the output images (not required)
		#	imgnames = raw_input("Image (thumbnail) names?: ")
		#	if len(imgnames) == 0: imgnames = galaxy

		if outfile == None: # Name of output script file
			outfile = input("Output (script) file name?: ")
			if len(outfile) == 0: outfile = galaxy+".sh" # if no name is given, just use the galaxy name
			if len(outfile.split(".")) < 2: outfile = outfile+".sh" # if .sh is not included in the name, add it


		### If using unique scales for each source ###
		if unique_scale == True: 
			print("This part of the code is not complete yet.")
			### Insert way of reading in the file that contains all of the unique scalings per source and output a DS9 bash file that'll scale accordingly
			### Currently, this looks at the unique scalings file and finds all sources in there in the DataFrame. It might be more useful to look for the scaling for each DataFrame source read in. To do that, may want to convert the scalings to a DataFrame and perform as search on it. (4/6/20: STILL TO DO)
			scalings = None
			try: 
				try: scalings  = np.genfromtxt(galaxy+"_scalings.txt", dtype=str).tolist() 
				except: 
					if ".txt" in scales: 
						scalings = np.genfromtxt(scales, dtype=str).tolist()
					else: scalings = np.genfromtxt(raw_input("Unique scalings file?: "), dtype=None).tolist()
			except: print("\nUnique scalings file not found.")
			#try: 
			print("Writing " + outfile + "...")
			with open(outfile, 'w') as f: 
				f.write("#! /bin/bash\necho \"Creating images of " + galaxy + " sources. Please wait...\"\n")
				# code for opening color image of galaxies
				f.write("ds9 -height "+str(imgsize[0])+" -width "+str(imgsize[1])+" -colorbar no \\\n-rgb -" + filterorder[0] + " -zscale " + colorfiles[0])	# green needs to be opened first, for alignment reasons
				f.write(" \\\n-" + filterorder[1] + " -zscale " + colorfiles[1])
				f.write(" \\\n-" + filterorder[2] + " -zscale " + colorfiles[2] + " \\\n") 
				# code for opening region files (if applicable)
				try: 
					for i in range(len(regions)): 
						f.write("-" + basefilter[i] + " -region load " + regions[i] + " ")
				except: pass;
				# For each of the scalings, look to see if the source is in the DataFrame
 				#  Making sure the scalings is in the right form
				if not isinstance(scalings[0], list): scalings = [scalings]
				# If so, save image. If not, pass.
				for k in range(len(scalings)):
					j = scalings[k]
					#try: 
					if isinstance(imgnames, list): tempname = str(imgnames[k])
					else: tempname = str(imgnames) + "_" + j[0] # if imgname is given, use it
					#except: tempname = j[0]

					# Search dataframe for the scaling ID
					temp = FindUnique(frame, idheader + " = " + j[0])
					
					# If ID is found, print pan script to .sh file
					if len(temp) > 0:
						f.write("\\\n-zoom to " + str(zoom) + " ")
						if coordsys == "img" or coordsys == "image":
							f.write("\\\n\\\n-" + base + " -pan to " + str(temp["x"].values[0]) + " " + str(temp["y"].values[0]) + " image " )
						elif coordsys == "galaxy" or coordsys == "fk5":
							f.write("\\\n-" + base + " -pan to " + str(temp["RA"].values[0]) + " " + str(temp["Dec"].values[0]) + " fk5 ")
						else:
							try: f.write("\\\n\\\n-" + base + " -pan to " + str(temp["x"].values[0]) + " " + str(temp["y"].values[0]) + " image " )
							except: f.write("\\\n-" + base + " -pan to " + str(temp["RA"].values[0]) + " " + str(temp["Dec"].values[0]) + " fk5 ")
						if j[1] == "-zscale" or j[1] == "zscale":
							f.write("\\\n-red -zscale -green -zscale -blue -zscale ")
						else: f.write("\\\n-red -scale limits 0 " + str(j[1]) + " -green -scale limits 0 " + str(j[2]) + " -blue -scale limits 0 " + str(j[3]) + " ")
						f.write("\\\n-saveimage " + filetype + " " + tempname +"."+filetype+" ")
						f.write("\\\n-zoom to " + str(env_zoom) + " ")
						f.write("\\\n-saveimage " + filetype + " " + tempname + "_env." + filetype + " ")
					else: pass;
				f.write("\\\n-exit\necho \"Done.\"")               
			f.close()
			print("DONE")
			#except: print("Error creating file.")



		### If not using unique scales for each source (default) ###
		else: 
			if scales == None:
				scales = []
				temp = "none" 
				while len(temp) > 0:
					temp = input("Color scales? (enter one at a time. Example: \"zscale\" or [\"red\", 10, 5, 2.5]. Press Enter when finished.)")
					if len(temp)>0: scales.append([x for x in re.split(",", temp)])
			
			if not isinstance(scales, list): scales = [scales] # Scales should be a list
			print("Writing " + outfile + "...")
			with open(outfile, 'w') as f: 
				f.write("#! /bin/bash\necho \"Creating images of " + galaxy + " sources. Please wait...\"\n")
				# code for opening color image of galaxies
				f.write("ds9 -height "+str(imgsize[0])+" -width "+str(imgsize[1])+" -colorbar no \\\n-rgb -" + filterorder[0] + " -zscale " + colorfiles[0])
				f.write(" \\\n-" + filterorder[1] + " -zscale " + colorfiles[1])
				f.write(" \\\n-" + filterorder[2] + " -zscale " + colorfiles[2] + " \\\n") 
				# code for opening region files
				try:
					for i in range(len(regions)): 
						f.write("-" + basefilter[i] + " -region load " + regions[i] + " ")
				except: pass;
				for i in range(len(frame)): 
					f.write("\\\n-zoom to " + str(zoom) + " ")
					if coordsys == "img" or coordsys == "image":
						f.write("\\\n-" + base + " -pan to " + str(frame["x"][i]) + " " + str(frame["y"][i]) + " image " )
					elif coordsys == "galaxy" or coordsys == "fk5":
						f.write("\\\n-" + base + " -pan to " + str(frame["RA"][i]) + " " + str(frame["Dec"][i]) + " fk5 ")
					else:
						try: f.write("\\\n-" + base + " -pan to " + str(frame["x"][i]) + " " + str(frame["y"][i]) + " image " )
						except: f.write("\\\n-" + base + " -pan to " + str(frame["RA"][i]) + " " + str(frame["Dec"][i]) + " fk5 ")
					for j in scales: 
						try: j = j.strip(",")
						except: j[0] = j[0].strip(",")
						if not isinstance(j, list): j = [j] # Convert to list

						# If a list of image names is given, use this. 
						# Else, if imgnames is a prefix, apply source numbers to image file name.
						# NOTE: SOMETHING IS GOING WRONG HERE WHEN UNIQUE SCALES AREN'T USED. WILL NEED TO LOOK INTO THIS LATER!
						if isinstance(imgnames, list): imgtemp = imgnames[i]
						else: imgtemp = imgnames+"%03i"%(i)

						if len(j) == 1:
							f.write("\\\n\\\n-red -zscale -green -zscale -blue -zscale ")
							pass;
						elif len(j) > 3: 
							f.write("\\\n\\\n-red -scale limits 0 "+str(j[1])+" -green -scale limits 0 "+str(j[2])+" -blue -scale limits 0 "+str(j[3])+" ")
							# In the case where multiple, non-unique scales are given, append current scale to imagenames
							if len(scales) > 1: imgtemp = [n + "_" + j[0] for n in imgtemp]
							else: pass;
							pass;
						f.write("\\\n-saveimage jpeg " + imgtemp + ".jpg ")
				f.write("\\\n-exit\necho \"Done.\"")               
			f.close()
			print("DONE") 
				
			


###-----------------------------------------------------------------------------------------------------


def WriteDS9_tile(frame=None, galaxy="M81", colorfiles=None, regions=None, scales="zscale", imgnames=None, imgsize=[305,305], outfile=None, unique_scale=False, coords=None, filetype="jpeg", zoom=8, env_zoom=2, coordsys=None, ids=None): 

	"""Arguments: (DataFrame, Galaxy name, Galaxy color .fits files in RGB order, Region files, Desired color scales, Output image name (number will be appended to the image), Output file name).

 Generates a bash script that will write a program for taking thumbnails/snapshots in DS9 of sources in a given DataFrame, but tiled by red, green, and blue! Will take in a list of region files for each frame, and different scales to use to adjust the RGB colors of the image. Scales should be input with the format ["scale name", redscale, greenscale, bluescale], or, if unique_scale is True, the name of the file with the unique scalings should be used. If zscale is being used, only ["zscale"] is needed. Parameter basefilter can be used to set which filter the region file should be aligned to, if not the default green (Ex. red used as the base filter in M81). The very first basefilter should be the main filter used for imaging, even if the region files are aligned to different filters. Filters should be called in the order they're intended to be used (default [red, green, blue]. If a different order is used, the order should be declared using the 'filterorder' parameter. UNDER CONSTRUCTION! """
	
	good = False

	# If only coordinates are given, build DataFrame from those.
	try: 
		if not frame: 
			if coordsys == "fk5" or coordsys == "wcs": 
				headers = ["ID", "RA", "Dec"]
			else: headers = ["ID", "x", "y"]
			if not ids: ids = np.arange(coords[0]).tolist()
			frame = BuildFrame(headers=headers, values=[ids, coords[0], coords[1]])
	except: pass;

	frame = frame.copy()

	#if colorfiles == None: # Name of the files for RGB frames needed
	#	while not good: 
	#		colorfiles = input("Names of .fits files for galaxy color map (red, green, blue; separate by commas): ")
	#		colorfiles = colorfiles.split(",")
	#		if len(colorfiles) == 3: good = True
	#		else: print("Not enough color files provided."); colorfiles = None; 

	if outfile == None: # Name of output script file
		outfile = input("Output (script) file name?: ")
		if len(outfile) == 0: outfile = galaxy+".sh" # if no name is given, just use the galaxy name
		if len(outfile.split(".")) < 2: outfile = outfile+".sh" # if .sh is not included in the name, add it


	print("Writing " + outfile + "...")
	with open(outfile, 'w') as f: 
		f.write("#! /bin/bash\necho \"Creating images of " + galaxy + " sources. Please wait...\"\n")
		# code for opening color image of galaxies
		f.write("ds9 -height "+str(imgsize[0])+" -width "+str(imgsize[1]*len(colorfiles))+" -colorbar no")

		# Open first frame and regions
		f.write(" \\\n-frame 1 " + colorfiles[0] + " -lock frame wcs -tile yes -tile column \\\n")
		try: 
			for i in range(len(regions[0])): 
				f.write("-region load " + regions[0][i] + " \\\n")
		except: pass;

		# Open the rest of the frames
		for i in range(1,len(colorfiles)):
			f.write("-frame " + str(i + 1) + " " + colorfiles[i] + " \\\n")
			try: 
				for j in range(len(regions[i])): 
					f.write("-region load " + regions[i][j] + " \\\n")
			except: pass;


		### If using unique scales for each source ###
		if unique_scale == True: 
			scalings = None
			try: 
				try: scalings  = np.genfromtxt(galaxy+"_scalings.txt", dtype=str) 
				except: 
					if ".txt" in scales: 
						scalings = np.genfromtxt(scales, dtype=str).tolist()
					else: scalings = np.genfromtxt(raw_input("Unique scalings file?: "), dtype=None).tolist()
			except: print("\nUnique scalings file not found.")
			
			# The code below assumes more than one unique scale is in the list
			# If it isn't, set scalings as a list within a list to make it work
			if not isinstance(scalings[0], list): scalings = [scalings]

			# For each of the scalings, look to see if the source is in the DataFrame
			# If so, save image. If not, pass.
			# In general, scalings is assumed to be a list of scalings with multiple sources included.
			# If only one source is in the unique scale file, will need to run another way
			for k in range(len(scalings)):
				j = scalings[k]
				if isinstance(imgnames, list): tempname = str(imgnames[k])
				else: tempname = str(imgnames) + "_" + j[0] # if imgname is given, use it

				# Search dataframe for the scaling ID
				temp = FindUnique(frame, "ID = " + j[0])

				# If ID is found, print pan script to .sh file
				if len(temp) > 0:
					f.write("\\\n-zoom to " + str(zoom) + " ")
					if coordsys == "img" or coordsys == "image":
						f.write("\\\n\\\n-frame 1 -pan to " + str(temp["x"].values[0]) + " " + str(temp["y"].values[0]) + " image " )
					elif coordsys == "galaxy" or coordsys == "fk5":
						f.write("\\\n-frame 1 -pan to " + str(temp["RA"].values[0]) + " " + str(temp["Dec"].values[0]) + " fk5 ")
					else:
						try: f.write("\\\n\\\n-frame 1 -pan to " + str(temp["x"].values[0]) + " " + str(temp["y"].values[0]) + " image " )
						except: f.write("\\\n-frame 1 -pan to " + str(temp["RA"].values[0]) + " " + str(temp["Dec"].values[0]) + " fk5 ")
					f.write("\\\n")
					if j[1] == "-zscale" or j[1] == "zscale":
						for i in range(1,len(colorfiles)+1): 
							f.write("-frame " + str(i) + " -zscale ")
					else: 
						for i in range(1,len(colorfiles)+1): 
							f.write("-frame " + str(i) + " -scale limits 0 " + str(j[i]) + " ")
					f.write("\\\n-saveimage " + filetype + " " + tempname +"."+filetype+" ")
					f.write("\\\n-zoom to " + str(env_zoom) + " ")
					f.write("\\\n-saveimage " + filetype + " " + tempname + "_env." + filetype + " ")
				else: pass;


		### If not using unique scales for each source (default) ###
		### NOTE: This part still needs to be edited to allow any number of colorfiles to be tiled, not just 3. (2/9/22) ###
		else: 
			if scales == None:
				scales = []
				temp = "none" 
				while len(temp) > 0:
					temp = input("Color scales? (enter one at a time. Example: \"zscale\" or [\"red\", 10, 5, 2.5]. Press Enter when finished.)")
					if len(temp)>0: scales.append([x for x in re.split(",", temp)])
			
			if not isinstance(scales, list): scales = [scales] # Scales should be a list
			for i in range(len(frame)): 
				f.write("\\\n-zoom to " + str(zoom) + " ")
				if coordsys == "img" or coordsys == "image":
					f.write("\\\n-frame 1 -pan to " + str(frame["x"][i]) + " " + str(frame["y"][i]) + " image " )
				elif coordsys == "galaxy" or coordsys == "fk5":
					f.write("\\\n-frame 1 -pan to " + str(frame["RA"][i]) + " " + str(frame["Dec"][i]) + " fk5 ")
				else:
					try: f.write("\\\n-frame 1 -pan to " + str(frame["x"][i]) + " " + str(frame["y"][i]) + " image " )
					except: f.write("\\\n-frame 1 -pan to " + str(frame["RA"][i]) + " " + str(frame["Dec"][i]) + " fk5 ")
				for j in scales: 
					try: j = j.strip(",")
					except: j[0] = j[0].strip(",")
					if not isinstance(j, list): j = [j] # Convert to list

					# If a list of image names is given, use this. 
					# Else, if imgnames is a prefix, apply source numbers to image file name.
					if len(imgnames) > 1: imgtemp = imgnames[i]
					else: imgtemp = imgnames+"%03i"%(i)

					if len(j) == 1:
						f.write("\\\n\\\n-frame 1 -zscale -frame 2 -zscale -frame 3 -zscale ")
						pass;
					elif len(j) > 3: 
						f.write("\\\n\\\n-frame 1 -scale limits 0 "+str(j[1])+" -frame 2 -scale limits 0 "+str(j[2])+" -frame 3 -scale limits 0 "+str(j[3])+" ")
						# In the case where multiple, non-unique scales are given, append current scale to imagenames
						if len(scales) > 1: imgtemp = [n + "_" + j[0] for n in imgtemp]
						else: pass;
						pass;
					f.write("\\\n-saveimage jpeg " + imgtemp + ".jpg ")

		### END WRITE TO FILE ###

		f.write("\\\n-exit\necho \"Done.\"")               
		f.close()
		print("DONE")
			
			


###-----------------------------------------------------------------------------------------------------

### UNDER CONSTRUCTION ###
# Need to figure out how I want to implement different variations on region files. 
def WriteReg(sources, outfile, coordsys=False, coordnames=False, idname=False, props=None, label=False, color="#FFC107", radius=3, radunit="pixel", showlabel=False, width=1, fontsize=10, bold=False, addshift=[0,0], savecoords=None, marker=None): 
	
	"""
	Writes a region file for all sources within the given DataFrame. Input (DataFrame, desired filename, coordinate system (image, fk5, etc), 
	region properties (region radius [including \" for arcsec], text, etc.) Arguments xcoord and ycoord refer to the header name of the 
	coordinates in the DataFrame.
	
	PARAMETERS
	----------
	sources     [DataFrame or list] :   Sources for which to plot the regions. Can be provided as either a 
	                                    DataFrame containing the headers [x,y] or [RA, Dec], or a list of 
	                                    coordinates in [[xcoords], [ycoords]] format. User can define the header 
	                                    of the coordinates with xcoord and ycoord.
	outfile     [str]               :   Name of the file to save to.
	coordsys    [str]               :   Defines the coordinate system to use in DS9 (e.g. image, fk5, etc.).
	coordnames  [list]              :   Name of the headers containing the x and y coordinates of each source, 
	                                    read in as a list in [xcoordname, ycoordname] format.
	idname      [str]               :   Name of the header containing the ID of each source. By default, checks 
	                                    whether the DataFrame contains a header called 'ID'. If not, it's assumed
	                                    no IDs are given for each source. 
	
	
	"""

	size=width
	
	if width:
		if not isinstance(width, list): width = "width=" + str(width)
		else: width = ["width=" + str(w) for w in width]

	if bold: bold=" bold"
	else: bold = " normal" 
	
	if coordnames: 
		xcoord = coordnames[0]
		ycoord = coordnames[1]
        
	### THIS NEEDS TO BE CLEANED
	if marker: # If a marker is defined
		if label: showlabel=True
		
		# Different forms of input needs to be treated differently.
		# Checking whether sources is a DataFrame or a list of coordinates
		
		if isinstance(sources, pd.DataFrame):
			# Setting up the possible coordinate systems based on inputs, or lack thereof.		
			sources = sources.copy()
			
			# Attempting to define the coordinate system and coordinate names, if they are
			# Not well-defined by the user. This assumes the coordsys is eeither fk5 or image.
			if not coordsys and not coordnames: 
				if "RA" in sources.columns.tolist():
					xcoord = "RA" 
					ycoord = "Dec"
					coordsys = "fk5"
					try: 
						rad_temp = float(radius)
						#radius = str(rad_temp * pixtoarcs) + "\""
					except: pass;
				else: 
					xcoord = "x" 
					ycoord = "y" 
					coordsys = "image"	
			elif coordsys == "fk5" and not coordnames: 
				xcoord = "RA" 
				ycoord = "Dec"	
				try: 
					rad_temp = float(radius)
					#radius = str(rad_temp * pixtoarcs) + "\""
				except: pass;
			elif coordsys == "image" and not coordnames: 
				xcoord = "x" 
				ycoord = "y" 	

			x_coords = sources[xcoord].values
			y_coords = sources[ycoord].values
			
			# if the label is given, use them as the ids
			if label: ids = label
			elif not label and showlabel: # If not and if showlabel is set to true, determine a good label
			    if not idname: # if no idname is given, search for values
			        try: ids = sources["ID"].values
			        except: ids = np.arange(0, len(x_coords))   # If no header found, number sources
			    else:   # if idname is defined, search for this header
			        try: ids = sources[idname].values
			        except: # if this header is not found, do not label
			            print("No header called", idname, "found. Sources will not be ID'd.")
			            ids = ""
			    
		elif isinstance(sources, list): 
			if len(np.asarray(sources)) == 2: 
				x_coords = sources[0]
				y_coords = sources[1]
			else: 
				x_coords = np.array(sources).T[0]
				y_coords = np.array(sources).T[1]
		elif len(re.split("\.", sources)) > 1: # if the sources are given as a filename, use GetCoords
			x_coords, y_coords = GetCoords(infile=sources) # retrieves coords from the file

		if not coordsys: # if coordsys not given, use simple check to assign
			if max(x_coords) > 500 or max(y_coords) > 500: coordsys = "image"
			else: coordsys= "fk5"

		
		# Making sure outfile has a proper file extension (default = .reg)
		temp = re.split("\.", outfile)
		if len(temp) == 1: outfile = temp[0] + ".reg"

		if not isinstance(radius, list): radius = [radius]
		if not isinstance(radunit, list): radunit = [radunit]*len(radius)

		for i in range(len(radius)): 
			if "arcs" in radunit[i]: radius[i] = str(radius[i]) + "\""
			if "arcm" in radunit[i]: radius[i] = str(radius[i]) + "\'"
			else: pass;


		if len(radius) == 1: radius = radius*len(x_coords) 

		### Attempting to add the props argument. More thought needed.
		#label = []
		#if props: 
		#	temp = re.split(",| ", props)
		#	if len(temp) < 2: temp = re.split(",", props)
		#	for i in temp:

		### UNDER CONSTRUCTION ###
		# props should contain different properties of the region file, such as color? 


		print("Saving " + outfile)
		with open(outfile, 'w') as f: 
			f.write("# Region file format: DS9 version 4.1\nglobal color=" + str(color) + " dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"+str(coordsys)+"\n")
			for i in range(len(radius)):	# In case each source in original region file has more than one radius, here we want to convert to only one radius.  
				x_coords[i] = float(x_coords[i]) + addshift[0]
				y_coords[i] = float(y_coords[i]) + addshift[1]

				# In case widths aren't the same for all sources
				if isinstance(width, list): w = width[i]
				else: w = width

				if showlabel == True and str(x_coords[i]) != "nan" and str(y_coords[i]) != "nan": 
					try:
						#if frame["ID"][i] != frame["ID"][i-1]: 
						#f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", 0.2\") #text={"+str(ids[i])+"}\n")
						f.write("point("+str(x_coords[i])+", "+str(y_coords[i])+") #text={"+str(ids[i])+"} font=\"helvetica " + str(fontsize) + bold + "\" "+width+"\n")
						#else: pass;
					except: 
						#f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", 0.2\")\n")
						f.write("point("+str(x_coords[i])+", "+str(y_coords[i])+") # point="+marker+" "+str(size)+"\n")
				elif showlabel == False and str(x_coords[i]) != "nan" and str(y_coords[i]) != "nan": 
					f.write("point("+str(x_coords[i])+", "+str(y_coords[i])+") # point="+marker+" "+str(size)+"\n")

		if savecoords: 
			with open(savecoords, "w") as f: 
				np.savetxt(f, np.column_stack([x_coords, y_coords]))
			print("Saved " + savecoords)
	else: # If a marker is not defined
		if label: showlabel=True
		
		# Different forms of input needs to be treated differently.
		# Not all coordinates will be inserted as a DataFrame
		if isinstance(sources, pd.DataFrame):
			# Setting up the possible coordinate systems based on inputs, or lack thereof.
			if not coordsys: 
				if "RA" in sources.columns.tolist():
					xcoord = "RA" 
					ycoord = "Dec"
					coordsys = "fk5"
					try: 
						rad_temp = float(radius)
						radius = str(rad_temp * pixtoarcs) + "\""
					except: pass;
				else: 
					xcoord = "x" 
					ycoord = "y" 
					coordsys = "image"	
			elif coordsys == "fk5": 
				xcoord = "RA" 
				ycoord = "Dec"	
				try: 
					rad_temp = float(radius)
					radius = str(rad_temp * pixtoarcs) + "\""
				except: pass;
			elif coordsys == "image" and not coordnames: 
				xcoord = "x" 
				ycoord = "y" 	

			x_coords = sources[xcoord].values
			y_coords = sources[ycoord].values
			
			# if the label is given, use them as the ids
			if label: ids = label
			elif not label and showlabel: # If not and if showlabel is set to true, determine a good label
			    if not idname: # if no idname is given, search for values
			        try: ids = sources["ID"].values
			        except: ids = np.arange(0, len(x_coords))   # If no header found, number sources
			    else:   # if idname is defined, search for this header
			        try: ids = sources[idname].values
			        except: # if this header is not found, do not label
			            print("No header called", idname, "found. Sources will not be ID'd.")
			            ids = ""
			            
		elif isinstance(sources, list): 
			if len(np.asarray(sources)) == 2: 
				x_coords = sources[0]
				y_coords = sources[1]
			else: 
				x_coords = np.array(sources).T[0]
				y_coords = np.array(sources).T[1]
				
			if label: ids = label
			elif not label and showlabel: ids = np.arange(0, len(x_coords))
			
		elif len(re.split("\.", sources)) > 1: # if the sources are given as a filename, use GetCoords
			x_coords, y_coords = GetCoords(infile=sources) # retrieves coords from the file

		if not coordsys: # if coordsys not given, use simple check to assign
			if max(x_coords) > 500 or max(y_coords) > 500: coordsys = "image"
			else: coordsys= "fk5"

		###--- BEGIN SAVING OUTPUT FILE ---###
		if not outfile: outfile = input("Output filename?: ")

		# Making sure outfile has a proper file extension (default = .reg)
		temp = re.split("\.", outfile)
		if len(temp) == 1: outfile = temp[0] + ".reg"

		if not isinstance(radius, list): radius = [radius]
		if not isinstance(radunit, list): radunit = [radunit]*len(radius)

		for i in range(len(radius)): 
			if "arcs" in radunit[i]: radius[i] = str(radius[i]) + "\""
			if "arcm" in radunit[i]: radius[i] = str(radius[i]) + "\'"
			else: pass;


		if len(radius) == 1: radius = radius*len(x_coords) 

		### Attempting to add the props argument. More thought needed.
		#label = []
		#if props: 
		#	temp = re.split(",| ", props)
		#	if len(temp) < 2: temp = re.split(",", props)
		#	for i in temp:

		### UNDER CONSTRUCTION ###
		# props should contain different properties of the region file, such as color? 

		print("Saving " + outfile)
		with open(outfile, 'w') as f: 
			f.write("# Region file format: DS9 version 4.1\nglobal color=" + str(color) + " dashlist=8 3 width=1 font=\"helvetica 10 "+bold+" roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"+str(coordsys)+"\n")
			for i in range(len(radius)):	# In case each source in original region file has more than one radius, here we want to convert to only one radius.  
				x_coords[i] = float(x_coords[i]) + addshift[0]
				y_coords[i] = float(y_coords[i]) + addshift[1]

				# In case widths aren't the same for all sources
				if isinstance(width, list): w = width[i]
				else: w = width

				if showlabel == True and str(x_coords[i]) != "nan" and str(y_coords[i]) != "nan": 
					try:
						#if frame["ID"][i] != frame["ID"][i-1]: 
						#f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", 0.2\") #text={"+str(ids[i])+"}\n")
						f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", " + str(radius[i]) + ") #text={"+str(ids[i])+"} font=\"helvetica " + str(fontsize) + bold + "\" "+width+"\n")
						#else: paøss;
					except: 
						#f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", 0.2\")\n")
						f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", " + str(radius[i]) + ") # "+w+"\n")
				elif showlabel == False and str(x_coords[i]) != "nan" and str(y_coords[i]) != "nan": 
					f.write("circle("+str(x_coords[i])+", "+str(y_coords[i])+", " + str(radius[i]) + ") # "+w+"\n")

		if savecoords: 
			with open(savecoords, "w") as f: 
				np.savetxt(f, np.column_stack([x_coords, y_coords]))
			print("Saved " + savecoords)


###-----------------------------------------------------------------------------------------------------

def WriteFilledReg(sources=None, coords=None, outfile=None, filename=None, coordsys="fk5", xcoord=None, ycoord=None, props=None, label=None, color="#FFC107", radius=3.0, radunit="pixel", showlabel=False, width=None, inner_width=1, outer_width=1, fontsize=10, bold=False, addshift=[0,0], savecoords=None, marker="circle", fill="#FFC107", outline="#FFC107", vertices=None, angle=0, size=None, inner_size=None, outer_size=None): 
	
	"""Writes a region file of filled markers for all sources within the given DataFrame. Input (DataFrame, desired filename, coordinate system (image, fk5, etc), region properties (region radius [including \" for arcsec], text, etc.). This can probably be combined with WriteReg above, with enough effort."""


	# The default marker type is a circle with some given radius
	# Also possible choices: triangle, square, diamond, 
	# Need to come up with some mathematical expression for each type of object. 
	# Use width as the outline width

	# Square markers are called box is DS9
	# Marker may be either circle, box, diamond, cross, x, arrow, or boxcircle
	if  marker == "square": marker = "box"
	elif marker == "+": marker = "cross"

	if bold: bold=" bold"
	else: bold = " normal" 

	# If a singular width is given, overwide the inner and outer widths
	# Outer width is defaulted to twice the inner width
	if width: 
		inner_width = width
		outer_width = width*2
	# Setting up the proper coordinate system and
	# converting radius into the proper units

	if coordsys == "fk5": 
		xcoord = "RA"
		ycoord = "Dec"
		try: 
			rad_temp = float(radius)
			radius = str(rad_temp * pixtoarcs) + "\""
		except: pass;
	elif coordsys == "image" and xcoord==None: 
		xcoord = "x" 
		ycoord = "y" 	

	# Different forms of input needs to be treated differently.
	# Not all coordinates will be inserted as a DataFrame
	if sources is not None:
		if isinstance(sources, pd.DataFrame):
			# Reading in the Dataframe and pulling coordinates
			sources = sources.copy()
			if not coords: 
				x_coords = sources[xcoord].values
				y_coords = sources[ycoord].values
			# If IDs are in DataFrame, pull
			try: 
				ids = sources["ID"].values
				for i in len(range(ids)): ids[i] = "text={"+ids[i]+"} font=\"helvetica " + str(fontsize) + bold + "\" "
			except: ids = ""
		# If sources is a list of coordinates, set coords
		elif isinstance(sources, list): 
			if len(np.asarray(sources)) == 2: 
				x_coords = sources[0]
				y_coords = sources[1]
			else: 
				x_coords = np.array(sources).T[0]
				y_coords = np.array(sources).T[1]
		# If the sources is a filename, use GetCoords
		elif len(re.split("\.", sources)) > 1: 
			x_coords, y_coords = GetCoords(infile=sources) # retrieves coords from the file

	# if coordinates are given separately, read them in
	if coords: 
		x_coords = coords[0]
		y_coords = coords[1]


	# If a label is given, set show label to true and ids to label
	# NOTE: this will override the automatic labeling set by the IDs in the Dataframe
	if label: 
		showlabel=True
		ids = []
		for i in label: 
			ids.append("text={" + i + "} font=\"helvetica " + str(fontsize) + bold + "\" ") 

	# We don't want anything in ids if showlabel=False
	if not showlabel: 
		ids = [""]*len(x_coords)
	else: pass; 

	# Setting up widths and marker sizes
	if not isinstance(width, list): width = [width]*len(x_coords)
	
	# Setting up the sizes. Inner and outer sizes can be independently controlled. 
	# If not given as size or inner/outer, check if given as radius, set that size.
	if not size: 
		if not inner_size and not outer_size: 
			try: 
				inner_size=int(radius)
				outer_size=int(radius)
			except:
				inner_size = width
				outer_size = width
	elif size: 
		inner_size = size
		outer_size = size

	# Width may sometimes be a list, so we will set size to always be a list
	if not isinstance(inner_size, list):
		inner_size = [inner_size]*len(x_coords) 
		outer_size = [outer_size]*len(x_coords) 

	###--- BEGIN SAVING OUTPUT FILE ---###
	if not outfile: outfile = input("Output filename?: ")

	# Making sure outfile has a proper file extension (default = .reg)
	temp = re.split("\.", outfile)
	if len(temp) == 1: outfile = temp[0] + ".reg"

	if not isinstance(radius, list): radius = [radius]*len(x_coords) 
	if not isinstance(radunit, list): radunit = [radunit]*len(radius)

	# radunits tells the units of the radius. This helps set the formatting 
	# of the radius as we add it to the list
	for i in range(len(radius)): 
		if "arcs" in radunit[i]: radius[i] = str(radius[i]) + "\""
		if "arcm" in radunit[i]: radius[i] = str(radius[i]) + "\'"
		else: pass;


	### --- WRITING THE FILE --- ###
	print("Saving " + outfile)
	with open(outfile, 'w') as f: 
		f.write("# Region file format: DS9 version 4.1\nglobal color=" + str(color) + " dashlist=8 3 width=1 font=\"helvetica 10 "+bold+" roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"+str(coordsys)+"\n")

		for i in range(len(radius)):	# In case each source in original region file has more than one radius, here we want to convert to only one radius.  

			# If a shift is given, add it
			x_coords[i] = float(x_coords[i]) + addshift[0]
			y_coords[i] = float(y_coords[i]) + addshift[1]

			# In case widths aren't the same for all sources
			if isinstance(inner_width, list): 
				iw = inner_width[i]
				ow = outer_width[i]
			else: 
				iw = inner_width
				ow = outer_width
			f.write("point("+str(x_coords[i])+", "+str(y_coords[i])+") #point="+marker+" "+ str(inner_size[i]) + " width=" + str(iw) + " color="+fill+ " " +ids[i]+ " \n")
			f.write("point("+str(x_coords[i])+", "+str(y_coords[i])+") #point="+marker+" "+ str(outer_size[i]) + " width=" + str(ow) + " color="+outline+ " " +ids[i]+ " \n")

	if savecoords: 
		with open(savecoords, "w") as f: 
			np.savetxt(f, np.column_stack([x_coords, y_coords]))
		print("Saved " + savecoords)

###-----------------------------------------------------------------------------------------------------

def WriteDonor(frame, galaxy=None, colorfiles=None, regions=None, scales=None, unique_scale=True, outfile=None, sourcecoords=None, donorcoords=None): #, imgnames=None, outfile=None): 

	"""Writes the bash script for imaging each of the given donor stars. These thumbnails will be useful for verifying donor star for each source.\nArguments: ()"""

	if galaxy == None: galaxy = input("Name of galaxy: ") # Name of the galaxy needed

	good = False

	frame = frame.copy()

	# Coordinates for the sources and the donors should be different. 
	# If not given, just use the coords from the DataFrame.
	# Not ideal, but it works
	if isinstance(sourcecoords, pd.DataFrame): 
		tempcoords = sourcecoords.copy()
		sourcecoords = []
		try: 
			for i in range(len(frame)): 
				temp = FindUnique(tempcoords, "ID = " + frame[ID][i])
				sourcecoords.append([temp[X][0], temp[Y][0]])
		except: 
			for i in range(len(frame)): 
				temp = FindUnique(tempcoords, "ID = " + frame[ID][i])
				sourcecoords.append([temp[RA][0], temp[Dec][0]])
	elif sourcecoords == None: 
		sourcecoords = []
		try: 
			for i in range(len(frame)): sourcecoords.append([frame[X][i], frame[Y][i]])
		except: 
			for i in range(len(frame)): sourcecoords.append([frame[RA][i], frame[Dec][i]])
	
	if isinstance(donorcoords, pd.DataFrame): 
		tempcoords = donorcoords.copy()
		donorcoords = []
		try: 
			for i in range(len(tempcoords)): donorcoords.append([tempcoords[X][i], tempcoords[Y][i]])
		except: 
			for i in range(len(tempcoords)): donorcoords.append([tempcoords[RA][i], tempcoords[Dec][i]])
	elif donorcoords == None: 
		donorcoords = []
		try: 
			for i in range(len(frame)): donorcoords.append([frame[X][i], frame[Y][i]])
		except: 
			for i in range(len(frame)): donorcoords.append([frame[RA][i], frame[Dec][i]])

	if colorfiles == None: # Name of the files for RGB frames needed
		while not good: 
			colorfiles = input("Names of .fits files for galaxy color map (red, green, blue; separate by commas): ")
			colorfiles = colorfiles.split(",")
			if len(colorfiles) == 3: good = True
			else: print("Not enough color files provided."); colorfiles = None; 
	
	# To draw region circles: -regions command circle x y rad # color=[color] dash=1
	if regions == None: #I n case region files were forgotten, ask user for input. 
		regions = input("Region files? (separate by commas, or Enter if none):") 
		regions = regions.split(",")
	if not isinstance(regions, list): regions = [regions]  # region files should be lists

	#if imgnames == None: # Names of the output images (not required)
	#	imgnames = nput("Image (thumbnail) names?: ")
	#	if len(imgnames) == 0: imgnames = galaxy

	if outfile == None: # Name of output script file
		#outfile = input("Output (script) file name?: ")
		outfile = galay + "_donors.sh"
		if len(outfile) == 0: outfile = galaxy+".sh" # if no name is given, just use the galaxy name
		if len(outfile.split(".")) < 2: outfile = outfile+".sh" # if .sh is not included in the name, add it
	
	### If using unique scales for each source ###
	if unique_scale == True: 
		# Reading in scalings from file. This will be an array that contains source IDs and scales.
		try: 
			try: scalings  = np.genfromtxt(galaxy+"_scalings.txt", dtype=None) 
			except: scalings = np.genfromtxt(raw_input("Unique scalings file?: "), dtype=None)
		except: print("\nUnique scalings file not found.")
		try: 
			print("Writing " + outfile + "...")
			with open(outfile, 'w') as f: 
				f.write("#! /bin/bash\necho \"Creating images of " + galaxy + " donors. Please wait...\"\n")
				tempid = "None"			
				donortemp = 0	
				for k in range(len(frame)): # For each of the objects in the frame
					if tempid == frame[ID][k]: donortemp = donortemp + 1
					else: donortemp = 0
					# code for opening color image of galaxies
					f.write("ds9 -height 200 -width 200 -colorbar no \\\n-rgb -red -zscale " + colorfiles[0])
					f.write(" \\\n-green -zscale " + colorfiles[1])
					f.write(" \\\n-blue -zscale " + colorfiles[2] + " \\\n") 
					tempscales = np.where(scalings.T[0] == frame[ID][k])[0][0]  # source scale
					tempscales = [x for x in scalings[tempscales]]  # pulling proper scale
					if tempscales[1] == "-zscale" or tempscales[1] == "zscale":
						j = [tempscales[0], " -zscale", " -zscale", " -zscale"]
					else: j = [tempscales[0], " -scale limits 0 " + str(tempscales[1]), " -scale limits 0 " + str(tempscales[2]), " -scale limits 0 " + str(tempscales[3])]
					f.write("-red" + str(j[1])+ " \\\n-green" + str(j[2]) + " \\\n-blue" +  str(j[3]) + " ")
					# code for opening region files
					for i in regions: 
						f.write(" \\\n-region load " + i + " ")
					f.write(" \\\n-regions command \"circle " + str(donorcoords[k][0]) + " " + str(donorcoords[k][1]) + " 3 # color=green dash=1 width=2\" ")
					f.write("\\\n-zoom to 4 ")
					try:
						temp = frame[X] 
						f.write("\\\n-pan to " + str(sourcecoords[k][0]) + " " + str(sourcecoords[k][1]) + " image " )
					except: f.write("\\\n-pan to " + str(sourcecoords[k][0]) + " " + str(sourcecoords[k][1]) + " fk5 ")
					f.write("\\\n-saveimage jpeg " + frame[ID][k] + "_donor" + str(donortemp) + ".jpg ")
					f.write("\\\n-exit \\\n\\\n")
					tempid = frame[ID][k]
					f.write("#! /bin/bash\n")
				f.write("\\\necho \"Done.\"")   
			f.close()
			print("DONE") 
		except: print("Error creating file.") 



	### If not using unique scales for each source (default) ###
	#else: 
	#	if scales == None:
	#		scales = []
	#		temp = "none" 
	#		while len(temp) > 0:
	#			temp = input("Color scales? (enter one at a time. Example: \"zscale\" or [\"red\", 10, 5, 2.5]. Press Enter when finished.)")
	#			if len(temp)>0: scales.append([x for x in temp.split()])
	#	
	#	if not isinstance(scales, list): scales = [scales] # Scales should be a list
	#	print "Writing " + outfile + "..."
	#	with open(outfile, 'w') as f: 
	#		f.write("#! /bin/bash\necho \"Creating images of " + galaxy + " sources. Please wait...\"\n")
	#		# code for opening color image of galaxies
	#		f.write("ds9 -height 305 -width 305 -colorbar no \\\n-rgb -red -zscale " + colorfiles[0])
	#		f.write(" \\\n-green -zscale " + colorfiles[1])
	#		f.write(" \\\n-blue -zscale " + colorfiles[2] + " \\\n") 
	#		# code for opening region files
	#		for i in regions: 
	#			f.write("-region load " + i + " ")
	#		f.write("\\\n-zoom to 4 ")
	#		for i in range(len(frame)): 
	#			try: f.write("\\\n-pan to " + str(frame["x"][i]) + " " + str(frame["y"][i]) + " image " )
	#			except: f.write("\\\n-pan to " + str(frame["RA"][i]) + " " + str(frame["Dec"][i]) + " fk5 ")
	#			for j in scales: 
	#				try: j = j.strip(",")
	#				except: j[0] = j[0].strip(",")
	#				if not isinstance(j, list): j = [j]
	#				if len(j) == 1:
	#					f.write("\\\n\\\n-red -zscale -green -zscale -blue -zscale ")
	#					f.write("\\\n-saveimage png " + imgnames+"%03i"%(i)+"_zscale.png ")
	#					pass;
	#				elif len(j) > 3: 
	#					f.write("\\\n\\\n-red -scale limits 0 "+str(j[1])+" -green -scale limits 0 "+str(j[2])+" -blue -scale limits 0 "+str(j[3])+" ")
	#					f.write("\\\n-saveimage png "+imgnames+"%03i"%(i)+"_"+j[0]+".png ")
	#					pass;
	#		f.write("\\\n-exit\necho \"Done.\"")               
	#	f.close()
	#	print "DONE" 
				
			

#### EXAMPLE CODE TO COPY FOR BASH SCRIPT WRITING ####
#print "Writing Long_shifted.sh"
#with open("Long_shifted.sh", 'w') as f: 
#    f.write("#! /bin/bash\necho \"Creating images of Long X-ray sources. Please wait...\"\n")
#    f.write("ds9 -height 305 -width 305 -colorbar no \\\n-rgb -red -zscale m83_red.fits "+ \
#            "\\\n-green -zscale m83_green.fits \\\n-blue -zscale m83_blue.fits \\\n"+ \
#            "-region load Long_cleaned.reg -region load Long_shifted.reg " +\
#            "-region load M83_sig1_raw.reg "+\
#            "-region load Long_CSC.reg \\\n-zoom to 2 ")
#    for i in range(len(Long_CSC_unique_in)):
#        try: 
#            f.write("\\\n-pan to " + str(Long_CSC_unique_in[X][i]) + " " + \
#                    str(Long_CSC_unique_in[Y][i]) + " image \\\n-saveimage png " + \
#                    "LongCSC"+"%03i"%(i)+"_zscale.png ")
#            f.write("\\\n\\\n-red -scale limits 0 20 -green -scale limits 0 10 -blue -scale limits 0 5 ")
#            f.write("\\\n-saveimage png LongCSC"+"%03i"%(i) + "_dark.png ")
#            f.write("\\\n\\\n-red -scale limits 0 4 -green -scale limits 0 2 -blue -scale limits 0 1 ")
#            f.write("\\\n-saveimage png LongCSC"+"%03i"%(i) + "_bright.png ")
#            f.write("\\\n\\\n-red -zscale -green -zscale -blue -zscale ")
#        except: pass;
#    f.write("\\\n-exit\necho \"Done.\"")               
#f.close()

#print "DONE"

###---------------------------------------------

def WriteScalings(sources=None, outfile="scalings.txt", scalings=None, default_scalings=["zscale","zscale","zscale"], regions=None, savescalings="autoscalings.txt", coords_header=['X','Y'], idheader="ID"): 
	"""Script for automatically writing a default unique scaling for each source in the input DataFrame based on the image coordinates 
	and input square regions. It is assumed that the regions get brighter closer to the center of the image and that regions are read in 
	going from the outer regions moving inward. The regions are read in as rectangles in the form [[xmin, xmax], [ymin, ymax]] in image 
	coordinates. Regions will be given as [red, green, blue]. If the source is found outside of all of the regions given, assumed to be 
	in zscale. Input scalings and regions will automatically be saved in a file with the default name 'autoscalings.txt'. This file can 
	be read in place of scalings and regions in situations where the same regions can be for the galaxy in all cases. This same file can
	then be called in as the scalings parameter, in place of a list of scalings. 

	The best steps are the following: 

	(1) Manually create an autoscaling file with a list of [redscale, greenscale, bluescale, xmin, xmax, ymin, ymax]; 
	(2) Run WriteScalings using the new autoscaling file name name as the 'scalings' argument to create a file with a 
	    list of unique location-based scalings for each source; 
	(3) Run WriteDS9 using unique_scale=True and the new unique scaling file as the 'scales' argument. Use the resulting 
	    .sh file to get scaled images of each source. 

	"""

	sourcex = sources[coords_header[0]].values.tolist()
	sourcey = sources[coords_header[1]].values.tolist()
	
	if ".txt" in scalings: 
		temp = np.genfromtxt(scalings, dtype=str)
		scalings = []
		regions = []
		for i in temp: 
			scalings.append([i[0], i[1], i[2]])
			regions.append([[i[3], i[4]],[i[5], i[6]]])	
	else: 
		print ("Auto-saving scalings as " + savescalings)
		with open(savescalings, "w") as f: 
			f.write("# Scalings and region coordinates readable by WriteScalings\n")
			for i in range(len(scalings)): 
				temp = ""
				for j in range(3): temp = temp + str(scalings[i][j]) + " "
				for j in range(2): temp = temp + str(regions[i][j][0]) + " " + str(regions[i][j][1]) + " "
				f.write(temp + "\n")
		 
	with open(outfile, "w") as f: 
		f.write("### Last update: " + str(datetime.date.today()))
		f.write("\n### Keeps track of the best unique color scaling for each source")
		f.write("\n# [0] ID, [1] Red scale, [2] Green scale, [3] Blue scale\n")
		# Looking at each source independently
		for s in range(len(sources)):
			# Default scale is zscale. If the source falls within a given region, change accordingly 
			red = default_scalings[0]
			green = default_scalings[1]
			blue = default_scalings[2]
			# Each region (reg) will be a list filled with two lists, [xmin, xmax] and [ymin, ymax]
			for i in range(len(regions)):
				reg = regions[i]
				scale = scalings[i]
				xmin = reg[0][0]
				xmax = reg[0][1]
				ymin = reg[1][0]
				ymax = reg[1][1]
				if float(xmin) < sourcex[s] < float(xmax) and float(ymin) < sourcey[s] < float(ymax):
					red = str(scale[0])
					green = str(scale[1])
					blue = str(scale[2])
				# Each of these scalings will be overwritten if the source appears in a later (i.e. inner) region
			temp = str(sources[idheader][s]) + " "  + str(red) + " " + str(green) + " " + str(blue)
			f.write(temp)
			f.write("\n")
	f.close()
	print (outfile + " saved\nDONE.")
				


###-----------------------------------------------------------------------------------------------------

def WriteFig(images, outfile=None, dimensions=(8,5), folder=None): #, imgnames=None, outfile=None): 

	"""For writing LaTex code for figure containing thumbnails. Dimensions are given as rows x columns"""

	images = images.copy()

	imgpertable = dimensions[0]*dimensions[1]  	# Number of images per table
	numtables = int(math.ceil(float(len(images))/float(imgpertable))) 	# Number of tables needed
	tables = [""]*numtables 	# a list for holding all of the tables generated
	img = 0 	# keeping track of which image is being added to the table

	try: 
		for i in images: images[i] = folder+"/"+images[i]
	except: pass;

	print("Writing " + outfile + "...\n")

	for i in range(numtables):	# This needs to be repeated for every table
		tables[i] = tables[i] + "\\begin{table}[ht]\n\t\\centering\n\t\\begin{tabular}{c"
		for j in range(dimensions[1]): tables[i] = tables[i] + "c"
		tables[i] = tables[i] + "}\t\t"
		try: 
			for j in range(dimensions[0]): 
				tables[i] = tables[i] + "\n\t\t"
				for k in range(dimensions[1]): 
					tables[i] = tables[i] + " & \\includegraphics[width=" + str(0.75/dimensions[1])+"\\textwidth]{"+images[img]+"}"
					img = img + 1
				tables[i] = tables[i] + " \\\\"
		except: tables[i] = tables[i] + "\\\\"
		tables[i] = tables[i] + "\n\t\\end{tabular}\n\t\\caption{}\n\t\\label{}\n\\end{table}\n"
		
	with open(outfile, 'w') as f: 
		for table in tables: 
			f.write(table + "\n")

	f.close()
	print("Done")

###-----------------------------------------------------------------------------------------------------

def WriteTable(frame, outfile=None, headers=None, dimensions=None): 

	"""For writing LaTex code for tables containing data properties, etc."""

	frame = frame.copy()

	# If no dimensions are given, use full df.
	if dimensions == None: dimensions = [frame.shape[0], frame.shape[1]]

	# If no headers are given, use all headers in the DataFrame
	if headers == None: headers = list(frame.columns.values) 
	else: dimensions[1] = len(headers)	# else, use only the given headers.

	numtables = int(math.ceil(float(len(frame))/float(dimensions[0])))
	tables = [""]*numtables
	source = 0
	# If headers are specified, only include those headers. 
	
	print("Writing " + outfile + "...\n")

	for i in range(numtables):	# This needs to be repeated for every table
		#tables[i] = tables[i] + "\\pagebreak\n\\begin{sidewaysfigure}\n\t\\begin{tabular}{c|" 
		tables[i] = tables[i] + "\n\\begin{center}\n\\afterpage{\n\t\\begin{longtable}[ht!]{c|"
		for j in range(dimensions[1]-1): tables[i] = tables[i] + "c"  # adding centered columns 
		tables[i] = tables[i] + "}\n\t\\hline\n\t"
		for head in headers: tables[i] = tables[i] + str(head) + " & "
		tables[i] = tables[i] + "\\\\ "
		try: 
			for j in range(dimensions[0]): 
				tables[i] = tables[i] + "\n\t"
				for head in headers: 
					if head != headers[-1]:
						try: tables[i] = tables[i] + str(round(frame[head][source], 3)) + " & "
						except: tables[i] = tables[i] + str(frame[head][source]) + " & "
					else: 
						try: tables[i] = tables[i] + str(round(frame[head][source], 3))
						except: tables[i] = tables[i] + str(frame[head][source])
				source = source + 1
				tables[i] = tables[i] + " \\\\ "
		except: tables[i] = tables[i] + "\\\\ "
		#tables[i] = tables[i] + "\n\t\\end{tabular}\n\\caption{}\n\\label{}\n\end{sidewaysfigure}\n"
		tables[i] = tables[i] + "\n\t\\end{longtable}\n\\clearpage\n}\n\\end{center}"
		
	with open(outfile, 'w') as f: 
		for table in tables: 
			f.write(table + "\n")

	f.close()
	print("Done")	




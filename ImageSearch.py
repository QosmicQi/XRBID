###################################################################################
##########		For finding and plotting HST images		########### 
##########		Last update: Dec. 3, 2024 			########### 
##########		Update desc: Normalizing spacing, parameters	########### 
##########		  						########### 
##########		  						########### 
###################################################################################

import warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import glob

import pyvo as vo

from astropy.io import fits
import astropy.coordinates as coord
from astropy.units import Quantity
# For downloading files
from astropy.utils.data import download_file
from astropy.wcs import WCS

from IPython.display import Image as ipImage, display

# There are a number of relatively unimportant warnings that show up, so for now, suppress them:
warnings.filterwarnings("ignore", module="astropy.io.votable.*")
warnings.filterwarnings("ignore", module="pyvo.utils.xml.*")

#-------------------------------------------------------------

def FindProducts(table, headers, criteria, exclude=False):

    """
    Given a table of data products from a search of a service (HLA), searches table for products 
    that meet specified criteria and returns a new table with just those products.

    PARAMETERS
    ------------------------
    table     [table] : Table containing service search results from pyVO
    headers   [list]  : List of headers on which the search criteria should be applied
    criteria  [list]  : List of criteria to apply to each header.
    exclude   [bool]  : List of bools indicating whether the criteria is exclusionary, rather 
	                    than inclusive (for removing certain keywords under headers)

    RETURNS
    ------------------------
    new_table [table] : Table containing only the rows of the original table that meet the search criteria
    indices    [list] : A list containing the indices of good sources from the original table.
                        Used to pull the correct image URL from the original table for plotting. 

    """

    tab_mask = [i for i in range(len(table))] # Starting off, the mask contains the full range of the table

    # Header and Criteria are assumed to be lists in the code, allowing multiple criteria to be entered
    if not isinstance(headers, list): headers=[headers]
    if not isinstance(criteria, list): criteria=[criteria]
    if not exclude: exclude = [False for i in range(len(headers))]

    for i,head in enumerate(headers):
        temp_table = table[tab_mask]
        # Keeping only the parts of the mask that meet the criteria
        if not exclude[i]: tab_mask = [tab_mask[j] for j,h in enumerate(temp_table[head]) if criteria[i] in h] 
        # tab_mask gradually gets smaller as it searches for rows where are criteria are met. 
	# Then, it returns the new, shortened table.
        else: tab_mask = [tab_mask[j] for j,h in enumerate(temp_table[head]) if criteria[i] not in h]

    new_table = table[tab_mask]

    return new_table, tab_mask  
    
#-------------------------------------------------------------

def FindHST(coords, galaxy, savefile, source_id="", search_rad=5, overwrite=False, plotall=False, savesearch=False):
    
    """
    Finds HST images of a source from the HLA, given the coordinates of that source. Sources 
    may either be read in as a single pair of coordinates (in units degrees), or as a DataFrame 
    containing the coordinates. FindHST() will prompt the user to decide which images from a 
    list of potential matches to save to the output file, savefile. Searches first for ACS/WFC, 
    then for WFC3/UVIS.
    

    PARAMETERS
    --------------------
    coords    [list]            : Coordinates of the source of interest in degrees [RA, Dec].
                                  May also be read in as a DataFrame for multiple sources.
    galaxy     [str]            : Galaxy name, for saving information about each source.
    savefile   [str]            : Name of file to save search results to, for easy recall.
    source_id  [str]            : Name of the source, to add to save file for easy recall.
    search_rad [float] (5)      : Search radius around source, in arc seconds.
    overwrite  [bool] (False)   : If true, overwrites savefile if it already exists. 
                                  Otherwise, appends results of the search to the end of the file.
    plotall    [bool] (False)   : If true, automatically plots all images without requesting from user.
    savesearch [str] (False)	: Saves the initial HLA search under the input filename. This allows the user 
				  to search through all search results at once to choose which they want.

    RETURNS
    --------------------
    Saves a CSV file with the name given as savefile containing the source ID, galaxy name, source 
    coordinates, filter, PropID (ObsID), ExpTime, and .FITS file URL.
    """

    # First checking if the savefile is already an existing file. If it isn't, set the file up. 
    # Otherwise, the results will append to the end of the existing file.
    if len(glob.glob(savefile)) == 0: 
        with open(savefile, "w") as f:
            f.write("Galaxy,SourceID,RA,Dec,Detector,Filter,ObsID,ExpTime,Filename,URL\n")
        f.close()
    else: 
        #if "y" in input("Overwrite " + savefile + "? (Selecting 'no' will append new search to current file): ").lower():
        if overwrite:
            with open(savefile, "w") as f:
                f.write("Galaxy,SourceID,RA,Dec,Detector,Filter,ObsID,ExpTime,Filename,URL\n")
            f.close()
            
    # Reading in the CSV file so that I can keep an eye on what has already been done, 
    # since this process may happen over multiple different iterations. 
    tempdf = pd.read_csv(savefile)

    # If the source is a DataFrame, then go source-by-source. 
    # Otherwise, treat as a regular list
    coords_original = coords.copy()
    
    if isinstance(coords_original, pd.DataFrame): 
        print(len(coords_original), "sources found.")
        # Search for the coordinates. If they cannot be found, prompt for proper the header input. 
        if "RA" in coords_original.columns.tolist() and "Dec" in coords_original.columns.tolist(): 
            coord_x = coords_original["RA"].tolist()
            coord_y = coords_original["Dec"].tolist()
        else: 
            head_x, head_y = input("Coordinate headers? (separated by comma): ").split(",")
            coord_x = coords_original[head_x].tolist()
            coord_y = coords_original[head_y].tolist()
            
        # Making sure each source has an ID, if not given
        if len(source_id) == 0:
            head_id = input("ID header name? ")
            ids = coords_original[head_id].tolist()
        elif isinstance(source_id, str): # assumes all sources have the same ID 
            ids = [source_id]*len(coords_original)
            
        # Run search for each source in the DataFrame
        for i in range(len(coords_original)): 
            # Checking if source already exists in current savefile
            # If it does, ask user whether they want to skip this source
            print("Finding images for", ids[i], "("+str(i+1), "of", str(len(coords_original)) + ")")
            if len(tempdf[tempdf["SourceID"] == ids[i]]) > 0:
                if 'n' in input(ids[i]+" already saved to file. Skip? (Y/n) ").lower(): 
                    FindHST_search(coords=[coord_x[i], coord_y[i]], galaxy=galaxy, savefile=savefile, source_id=ids[i], search_rad=search_rad)
                else: pass;
            else: FindHST_search(coords=[coord_x[i], coord_y[i]], galaxy=galaxy, savefile=savefile, source_id=ids[i], search_rad=search_rad, plotall=plotall, savesearch=savesearch)
    else: 
        FindHST_search(coords_original, galaxy, savefile, source_id, search_rad, plotall=plotall, savesearch=savesearch)
        
#-------------------------------------------------------------

def FindHST_search(coords, galaxy, savefile, source_id="", search_rad=5, plotall=False, savesearch=False): 
    
    """ An auxiliary function for FindHST(), which actually performs the search and save. """
            
    # Setting up the search using HLA as the defined serivce
    hst_services = vo.regsearch(servicetype='image',
                                waveband='optical',
                                keywords=["hla"])
                                
    # Defining the coordinates as WCS sky coordinates
    source_coords = coord.SkyCoord(coords[0],coords[1],frame="fk5",unit="degree")

    print("Searching HLA...")
    # Results of an image search around source coordinates
    image_search = hst_services[0].search(pos=source_coords,
                                          size=Quantity(search_rad, unit="arcsecond"),
                                          format='image/fits')

    # Saving search results to a searchable table
    image_table = image_search.to_table()

    # Searching the HLA for good ACS images
    good_acs_table, acs_inds = FindProducts(image_table, headers=["Detector", "Spectral_Elt", "Spectral_Elt"],
               criteria=["ACS/WFC", "detection", "/"],
               exclude=[False, True, True])
    good_acs_table['index'] = [ind for ind in range(len(good_acs_table))] #adding indices to output

    # Saving the search results to a file, if prompted
    if savesearch:
        with open(savesearch, "a") as f:
          for ind in [ind for ind in range(len(good_acs_table))]: 
            select_good_table = good_acs_table[ind]
            select_ind = acs_inds[ind]              # Index of the ACS file in the original image table
            select_image = image_search[select_ind]
            select_url = select_image.getdataurl()  # url of the selected image
            select_str = ",".join((str(ind),galaxy,source_id,str(coords[0]),str(coords[1]),
                                  str(select_good_table["Detector"]), str(select_good_table["Spectral_Elt"]),
                                  str(select_good_table["PropID"]), str(select_good_table["ExpTime"]), 
                                  str(select_good_table["filename"]), str(select_url))) + "\n"
            f.write(select_str)
        print(savesearch, "saved.")
        f.close()

    if len(good_acs_table) > 0:
        print("\n\nMatching ACS/WFC Observations:")
        display(good_acs_table["index","filename", "PropID", "Detector", "Spectral_Elt","ExpTime"])

        # Prompts user if they want to plot each of the images listed above. If yes, plot.
        if plotall or "y" in input("Plot images? (y/N): ").lower():
            for i in range(len(good_acs_table)):
                print(i, good_acs_table["Title"][i])
                hdu = fits.open(image_search[acs_inds[i]].getdataurl())
                #hdu_list.info()
                try:
                    w = WCS(hdu[1].header)
                    plt.figure(figsize=(3,3))
                    ax = plt.subplot(projection=w)
                    plt.imshow(hdu[1].data, cmap='gray', origin='lower',
                              norm=matplotlib.colors.Normalize(vmin=-0.04, vmax=.5))
                    plt.scatter(coords[0], coords[1], \
                                    marker="x", s=50, lw=2, color="red",\
                                    transform=ax.get_transform('world'))
                    plt.show()
                except: print("Plotting failed.\n")
        ### End of plotting

        # Prompts user to select the images to save to a file
        # Saving the Galaxy, source coordinates, filter, PropID, ExpTime, filename, and url to a file
        good_select = input("Input index of FITS files you wish to save (separated by commas):").split(",")
        
        with open(savefile, "a") as f:
          for ind in good_select: 
            ind = int(ind.strip())                  # index of the good FITS files
            select_good_table = good_acs_table[ind]
            select_ind = acs_inds[ind]              # Index of the ACS file in the original image table
            select_image = image_search[select_ind]
            select_url = select_image.getdataurl()  # url of the selected image
            select_str = ",".join((galaxy,source_id,str(coords[0]),str(coords[1]),
                                  str(select_good_table["Detector"]), str(select_good_table["Spectral_Elt"]),
                                  str(select_good_table["PropID"]), str(select_good_table["ExpTime"]), 
                                  str(select_good_table["filename"]), str(select_url))) + "\n"
            f.write(select_str)
        f.close()

    else: print("\n\nNo ACS/WFC results.")


    # Searching the HLA for good WFC3 images
    good_wfc3_table, wfc3_inds = FindProducts(image_table, headers=["Detector", "Spectral_Elt", "Spectral_Elt"],
               criteria=["WFC3/UVIS", "detection", "/"],
               exclude=[False, True, True])
    good_wfc3_table['index'] = [ind for ind in range(len(good_wfc3_table))] #adding indices to output

    # Saving the search results to a file, if prompted
    if savesearch:
        with open(savesearch, "a") as f:
          for ind in [ind for ind in range(len(good_wfc3_table))]: 
            select_good_table = good_wfc3_table[ind]
            select_ind = wfc3_inds[ind]              # Index of the WFC3 file in the original image table
            select_image = image_search[select_ind]
            select_url = select_image.getdataurl()  # url of the selected image
            select_str = ",".join((str(ind),galaxy,source_id,str(coords[0]),str(coords[1]),
                                  str(select_good_table["Detector"]), str(select_good_table["Spectral_Elt"]),
                                  str(select_good_table["PropID"]), str(select_good_table["ExpTime"]), 
                                  str(select_good_table["filename"]), str(select_url))) + "\n"
            f.write(select_str)
        print(savesearch, "updated.")
        f.close()
    
    if len(good_wfc3_table) > 0:
        print("\n\nMatching WFC3/UVIS Observations:")
        display(good_wfc3_table["index","filename", "PropID", "Detector", "Spectral_Elt","ExpTime"])

        # Prompts user if they want to plot each of the images listed above. If yes, plot.
        if plotall or "y" in input("Plot images? ").lower():
            for i in range(len(good_wfc3_table)):
                #print(i, good_wfc3_table["Title"][i])
                print(i, good_wfc3_table["filename"][i])
                hdu = fits.open(image_search[wfc3_inds[i]].getdataurl())
                #hdu_list.info()
                try:
                    w = WCS(hdu[1].header)
                    plt.figure(figsize=(3,3))
                    ax = plt.subplot(projection=w)
                    plt.imshow(hdu[1].data, cmap='gray', origin='lower',
                              norm=matplotlib.colors.Normalize(vmin=-0.04, vmax=.5))
                    plt.scatter(coords[0], coords[1], \
                                    marker="x", s=50, lw=2, color="red",\
                                    transform=ax.get_transform('world'))
                    plt.show()
                except: print("Plotting failed.\n")
        ### End of plotting

        # Prompts user to select the images to save to a file
        # Saving the Galaxy, source coordinates, filter, PropID, ExpTime, filename, and url to a file
        good_select = input("Input index of FITS files you wish to save (separated by commas):").split(",")
        
        with open(savefile, "a") as f:
            for ind in good_select: 
                ind = int(ind.strip())                  # index of the good FITS files
                select_good_table = good_wfc3_table[ind]
                select_ind = wfc3_inds[ind]              # Index of the ACS file in the original image table
                select_image = image_search[select_ind]
                select_url = select_image.getdataurl()  # url of the selected image
                select_str = ",".join((galaxy,source_id,str(coords[0]),str(coords[1]),
                                      str(select_good_table["Detector"]), str(select_good_table["Spectral_Elt"]),
                                      str(select_good_table["PropID"]), str(select_good_table["ExpTime"]), 
                                      str(select_good_table["filename"]), str(select_url))) + "\n"
                f.write(select_str)
        f.close()   
    else: print("\n\nNo WFC3/UVIS results.")


#-------------------------------------------------------------

def GetFITS(infile, source_id, detector, filter): 
    
    """
    Calls in the appropriate FITS file for the desired source, detector
    
    PARAMETERS
    ------------------------
    infile      [str] : Name of the file containing the observation log, as written by FindHST(). 
    source_id   [str] : Given ID of the source of interest. 
    detector    [str] : Name of the detector (ACS/WFC or WFC3/UVIS). 
    filter      [str] : Name of the filter of interest (e.g. F814W, F606W, etc.). 
    
    RETURNS
    ------------------------
    Returns the HDU containing the FITS file of the image within which the source falls. 

	"""

    filename = infile
    Obslog = pd.read_csv(filename)

    # Finding matches in the observation log
    select_obs = Obslog[Obslog["SourceID"] == source_id].reset_index(drop=True)
    select_obs = select_obs[select_obs["Detector"] == detector].reset_index(drop=True)
    select_obs = select_obs[select_obs["Filter"] == filter].reset_index(drop=True)
    
    # If more than one match exists, prompt user to select the appropriate file
    if len(select_obs) > 1: 
        display(select_obs)
        select_ind = int(input("Input index of the desired observation:"))
    else: select_ind = 0

    hdu = fits.open(select_obs["URL"][select_ind])

    return hdu
    
###-----------------------------------------------------------------------------------------------------

def MakePostage(gal, df, fits_dir, main_dir, reg_dir, instrument=None, id_header="ID", autolabel=False, obslog=False, img_width = 5., img_format="jpg", annotation=False, Clusters=None):

	"""
	NEEDS TO BE UPDATED TO GENERALIZE THE CODE (currently ULX-specific language used)
	Generates HST postage-stamp images of X-ray sources, to be used in a larger mosaic.
	These images are, by default, 5 arcseconds in width. This code checks to see if there is an Obslog
	file, which is produced by FindHST(). If there is no Obslog file, looks for the FITS files for 
	the galaxy of interest within the fits directory. 

 	PARAMETERS
	----------
	gal		[str]           :   Galaxy name.
	df		[pd.DataFrame]	:   A pandas DataFrame filled with the X-ray sources of interest.
	fits_dir	[str]           :   Directory of the FITS files
	main_dir	[str]           :   Main directory, to which the images will be saved
	reg_dir		[str]           :   Directory of the region files
	instrument	[str] (None)    :   Name of the detector of interest. If none is given, find in the obslog, 
					    fits file name, or prompt user input.
	id_header	[str] ('ID')    :   The name of the header in which the X-ray source IDs may be found in the 
					    source DataFrame.
	autolabel	[bool] (False)  :   If True, labels the image by the ID of the source. Otherwise, requests 
					    user input.
	obslog		[bool] (False)  :   Declares whether the galaxy has an observation log (relevant for newer
                                  	    galaxies added to analysis). If False, searches directory for relevant 
					    fits files
	img_width	[float] (5.)    :   Image width in arcseconds.
	img_format	[str] ('jpeg')  :   Format of the image file to save.
	annotation	[str] (False)   :   Annotation to add to the image (i.e. galaxy name).
	Clusters	[pd.DataFrame]	:   Dataframe containing clusters to highlight.

	RETURNS
	---------
	Saves an image file with the namingi convention ULX_[ULX ID]_[INSTRUMENT SHORTHAND].[image format]

	"""

	sources = df.copy()
	Temp_sources = Find(sources, "Galaxy = " + gal)
	tempulxs = FindUnique(Temp_sources, header=id_header)[id_header].tolist()

	# Plotting each ULX in said galaxy
	for u,tempulx in enumerate(tempulxs):

		print("Working on",tempulx)
		# If there is an observation log, read necessary information from that log. 
		#Ask user to input RGB filters from the list that appears.

		if obslog:
			cd(main_dir)
			Temp_obs = pd.read_csv("Obslog_"+gal+".csv")
			Temp_ulx_obs = Find(Temp_obs, "SourceID = " + tempulx) # pull the appropriate observation file

			# Saving and sorting list of filter names. If there are only 2 filters, add None to list
			display(Temp_ulx_obs)

			time.sleep(0.5) # pause before requesting input

			# Requesting user to choose
			fits_inds = input("Indices of correct RGB filters (comma, no space). If no filter available, input as 						  'None':").split(",")
			fits_inds = [int(value) if value != "None" else value for value in fits_inds]

			# Compiling filter names
			tempi, tempv, tempb = [Temp_ulx_obs["Filter"][value] if value!="None" else "—" for value in fits_inds]
			print(tempi, tempv, tempb)

			cd(fits_dir)
			# Setting default header for WCS info
			# By default, use V-band. Else, use I-band
			if tempv != "-":
				hdu_head = fits.open(Temp_ulx_obs["URL"][fits_inds[1]])
				goodind = 1
			else:
				hdu_head = fits.open(Temp_ulx_obs["URL"][fits_inds[0]])
				goodind = 0

			# Opening fits files. If none is given, set to an array of zeros.
			print("Reading in FITS files...")
			if tempi != "—":
				data_r = fits.getdata(Temp_ulx_obs["URL"][fits_inds[0]])
				if data_r.shape != hdu_head[1].data.shape: # Making sure all images are the same shape
					print("Reshaping I-band image...")
					data_r, _ = reproject_interp(fits.open(Temp_ulx_obs["URL"][fits_inds[0]])[1], \
								     hdu_head[1].header)
			else: data_r = np.zeros_like(hdu_head[1].data)

			if tempv != "—": # V-band will always be the reference filter if available.
				data_g = fits.getdata(Temp_ulx_obs["URL"][fits_inds[1]])
			else: data_g = np.zeros_like(hdu_head[1].data)

			if tempb != "—":
				data_b = fits.getdata(Temp_ulx_obs["URL"][fits_inds[2]])
				if data_b.shape != hdu_head[1].data.shape: # Making sure all images are the same shape
					print("Reshaping B-band image...")
					data_b, _ = reproject_interp(fits.open(Temp_ulx_obs["URL"][fits_inds[2]])[1], \
								     hdu_head[1].header)
			else: data_b = np.zeros_like(hdu_head[1].data)

			# set up the instrument name
			if not instrument: instrument = Temp_ulx_obs["Detector"][fits_inds[goodind]]
			else:
				if "acs" in instrument.lower(): instrument="ACS/WFC"
				else: instrument = "WFC3/UVIS"

			cd(main_dir)
			##### End first part of IF statement #####

		else: # If there is no obslog, search for FITS files directly
			cd(fits_dir)
			fits_list = glob.glob("*"+gal+"*.fits") # pulling file names
			fits_list.sort() # sorting

			# Printing FITS files found related to galaxy
			for ind,f in enumerate(fits_list): print(ind,f)

			time.sleep(0.5) # pause before requesting input

			# Requesting user to choose
			fits_inds = input("Indices of correct RGB filters (comma, no space). If no filter available, input as 'None':").split(",")
			fits_inds = [int(value) if value != "None" else value for value in fits_inds]
			# Compiling filter names
			tempi, tempv, tempb = [fits_list[int(value)].split("_")[1].split(".")[0].upper() if value!="None" \
					       else "—" for value in fits_inds]
			print(tempi, tempv, tempb)

			# Setting default header for WCS info
			# By default, use V-band. Else, use I-band
			if tempv != "—":
				hdu_head = fits.open(fits_list[fits_inds[1]])
				goodind = 1 # index of the header to use for the base hdu
			else:
				hdu_head = fits.open(fits_list[fits_inds[0]])
				goodind = 0 # index of the header to use for the base hdu

			# Opening fits files. If none is given, set to zero.
			cd(fits_dir)
			print("Reading in FITS files...")
			if tempi != "—":
				data_r = fits.getdata(fits_list[fits_inds[0]])
				if data_r.shape != hdu_head[1].data.shape: # Making sure all images are the same shape
					print("Reshaping I-band image...")
					data_r, _ = reproject_interp(fits.open(fits_list[fits_inds[0]])[1], hdu_head[1].header)
			else: data_r = np.zeros_like(hdu_head[1].data)

			if tempv != "—": data_g = fits.getdata(fits_list[fits_inds[1]])
			else: data_g = np.zeros_like(hdu_head[1].data)

			if tempb != "—":
				data_b = fits.getdata(fits_list[fits_inds[2]])
				if data_b.shape != hdu_head[1].data.shape: # Making sure all images are the same shape
					print("Reshaping B-band image...")
					data_b, _ = reproject_interp(fits.open(fits_list[fits_inds[1]])[1], hdu_head[1].header)
			else: data_b = np.zeros_like(hdu_head[1].data)

			# set up the instrument name
			if not instrument:
				if "acs" not in Temp_ulx_obs["URL"][fits_inds[goodind]].lower() and "wfc3" not in \
				Temp_ulx_obs["URL"][fits_inds[goodind]].lower(): instrument = input("Instrument/Detector?")

			# Standardizing the instrument format
			if "acs" in instrument.lower(): instrument="ACS/WFC"
			else: instrument = "WFC3/UVIS"
			cd(main_dir)

		##### END IF STATEMENT #####


		# Once the fits files are all set up, can start plotting.
		# First will plot the full HST FOV and ask for a color adjustment
		# Then, will post a close-up of the source
		cd(main_dir)
		Temp_ulx = Find(Temp_sources, "ID = " + tempulx)
		tempra = Temp_ulx["RA"][0]
		tempdec = Temp_ulx["Dec"][0]

		# Coordinates of the source
		wcs = WCS(hdu_head[1].header)
		tempra_pix,tempdec_pix = wcs.world_to_pixel_values(tempra, tempdec)

		adjustments = [1,1,1] # Factor adjustment of RGB
		m = 0.1 # the minimum mag cutoff of the RGB image
		M = 1   # The maximum mag cuttoff of the RGB image

		# Putting together the color image
		# No longer using Lupton RGB; instead, used my homemade log scaler
		rgb = scaled_log_intensity(data_r, data_g, data_b,  m=m, M=M)

		# First plot full galaxy to get the colors right
		plt.figure(figsize=(4,4))
		ax = plt.subplot(projection=wcs) # Setting the plot to the coordinate system
		plt.imshow(rgb)
		plt.axis("off")
		plt.show()
		plt.pause(1)

		while "y" in input("Adjust colors? ").lower():
			# If the colors are wrong, input factor by which to adjust brightness of each filter as fraction of current values
			time.sleep(0.1)
			adjustments = input("Input factor of adjustment (r,g,b): ").split(",")
			new_r = data_r*float(adjustments[0])
			new_g = data_g*float(adjustments[1])
			new_b = data_b*float(adjustments[2])
			rgb = scaled_log_intensity(new_r, new_g, new_b, m=m, M=M)
			plt.figure(figsize=(4,4))
			ax = plt.subplot(projection=wcs) # Setting the plot to the coordinate system

			plt.imshow(rgb)
			plt.axis("off")
			plt.show()
			plt.pause(1)

		# Updating filters
		new_r = data_r*float(adjustments[0])
		new_g = data_g*float(adjustments[1])
		new_b = data_b*float(adjustments[2])

		print("Plotting source...")
		rgb = scaled_log_intensity(new_r, new_g, new_b, m=m, M=M)
		plt.figure(figsize=(4,4))
		ax = plt.subplot(projection=wcs) # Setting the plot to the coordinate system
		plt.imshow(rgb)

		ax.set_xlim(tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
			    tempra_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))
		ax.set_ylim(tempdec_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
			    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))

		plt.axis("off")
		plt.show()
		plt.pause(1)

		print("Current min and max brightness:",m,M)
		while "y" in input("Adjust brightness?").lower():

			m,M = input("Input minimum and maximum brightness (m,M): ").split(",")
			m = float(m)
			M = float(M)
			rgb = scaled_log_intensity(new_r, new_g, new_b, m=m, M=M)

			plt.figure(figsize=(4,4))
			ax = plt.subplot(projection=wcs) # Setting the plot to the coordinate system
			plt.imshow(rgb)
			ax.set_xlim(tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
				    tempra_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))
			ax.set_ylim(tempdec_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
				    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))
			plt.axis("off")
			plt.show()
			plt.pause(0.5)

		# One last color adjustment
		while "y" in input("Adjust colors? ").lower():
			# If the colors are wrong, input factor by which to adjust brightness of each filter
			adjustments = input("Input factor of adjustment (r,g,b): ").split(",")

			new_r = data_r*float(adjustments[0])
			new_g = data_g*float(adjustments[1])
			new_b = data_b*float(adjustments[2])
			rgb = scaled_log_intensity(new_r, new_g, new_b, m=m, M=M)

			plt.figure(figsize=(4,4))
			ax = plt.subplot(projection=wcs) # Setting the plot to the coordinate system
			plt.imshow(rgb)
			ax.set_xlim(tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
				    tempra_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))
			ax.set_ylim(tempdec_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
				    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))
			plt.axis("off")
			plt.show()
			plt.pause(0.5)


		# Plotting ULX one last time and saving
		new_r = data_r*float(adjustments[0])
		new_g = data_g*float(adjustments[1])
		new_b = data_b*float(adjustments[2])
		rgb = scaled_log_intensity(new_r, new_g, new_b, m=m, M=M)

		plt.figure(figsize=(4,4))
		ax = plt.subplot(projection=wcs) # Setting the plot to the coordinate system
		plt.imshow(rgb)

		# Plotting filters used in image
		ax.annotate(instrument, (tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"])+1,
			    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"])-7),
			    color="white", fontsize=12, weight="bold",
			    path_effects=[pe.withStroke(linewidth=2, foreground="black")])

		ax.annotate(tempi, (tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"])+1,
			    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"])-15),
			    color="red", fontsize=12, weight="bold",
			    path_effects=[pe.withStroke(linewidth=2, foreground="white")])

		ax.annotate(tempv, (tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"])+1,
			    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"])-23),
			    color="limegreen", fontsize=12, weight="bold",
			    path_effects=[pe.withStroke(linewidth=2, foreground="black")])

		ax.annotate(tempb, (tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"])+1,
			    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"])-31),
			    color="blue", fontsize=12, weight="bold",
			    path_effects=[pe.withStroke(linewidth=2, foreground="white")])

		# Adding additional annotation
		if annotation:
			ax.annotate(annotation, (tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"])+2,
				    tempdec_pix-img_width/2./float(hdu_head[0].header["D001SCAL"])+2),
				    color="white", fontsize=15, weight="bold",
				    path_effects=[pe.withStroke(linewidth=2, foreground="black")])


		ulx_1sig = CircularAperture(np.array([tempra_pix, tempdec_pix]), r=Temp_ulx["1sig"][0]/float(hdu_head[0].header["D001SCAL"]))
		ulx_1sig.plot(color="red", lw=2)
		ulx_2sig = CircularAperture(np.array([tempra_pix, tempdec_pix]), r=Temp_ulx["2sig"][0]/float(hdu_head[0].header["D001SCAL"]))
		ulx_2sig.plot(color="red", lw=3)

		# THE DIMENSIONS OF EACH IMAGE ARE 5 arcseconds by 5 arcseconds #
		ax.set_xlim(tempra_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
			    tempra_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))
		ax.set_ylim(tempdec_pix-img_width/2./float(hdu_head[0].header["D001SCAL"]),
			    tempdec_pix+img_width/2./float(hdu_head[0].header["D001SCAL"]))

		cd(reg_dir)
		# Reading and plotting in region files for candidate stars
		reg_files = glob.glob("*"+tempulx+"*_stars*.reg")

		if len(reg_files) > 1:
			for ind,regf in enumerate(reg_files): print(ind,regf)
			reg_ind = int(input("Index of correct region file: "))
			reg_file = reg_files[reg_ind]
		elif len(reg_files) == 1: reg_file = reg_files[0]

		if len(reg_files) > 0:
			starcoords = GetCoords(reg_file, verbose=False)
			starids = GetIDs(reg_file, verbose=False)

			cd(main_dir)
			# Reading in and plotting candidate clusters
			if len(Clusters) > 0:
				tempclust = Find(Clusters, ["ULX = " + tempulx, "Instrument = " + instrument.split("/")[0].upper()])
				tempclust_ids = tempclust["Star"].tolist()
				tempclust_ids = [int(clust) for clust in tempclust_ids]
			else: tempclust_ids = [None]

			# Checking if each source is a cluster
			for ind,s in enumerate(starids):
				if int(s) in tempclust_ids: # if the source is a cluster, plot as a dashed circle
					plt.scatter(starcoords[0][ind], starcoords[1][ind], facecolors='none',
						    edgecolors="red", linestyle="--", s=450,lw=2)
				else: # otherwise, use a solid circle
					plt.scatter(starcoords[0][ind], starcoords[1][ind], facecolors='none', edgecolors="red", s=250,lw=1)

				ax.annotate(s, (starcoords[0][ind]-2.5, starcoords[1][ind]+8), color="white", fontsize=12, weight="bold",
					    path_effects=[pe.withStroke(linewidth=1.5, foreground="black")])

		else: pass;

		cd(main_dir)
		plt.axis("off")
		if autolabel == True: plt.title(tempulx,fontsize=20)
		elif autolabel != False: plt.title(autolabel, fontsize=20)
		else:
			print(tempulx)
			plt.title(input("Label? "), fontsize=20)

		plt.savefig("ULX_" + tempulx + "_"+instrument.split("/")[0].lower() + "." + img_format,
			    bbox_inches="tight", dpi=300)
		plt.show()
		plt.pause(0.5)

		print("ULX_" + tempulx + "_"+instrument.split("/")[0].lower() + "." + img_format, "saved")

#-------------------------------------------------------------

def PlotField(infile, fields=False, source_id=False): 
    
    """
    Creates a color image from the FITS files associated with a specific HST field of view. 
    This can be used to plot the entire field of view, instead of focusing in on a single
    source as MakePostage() does. 
    
    PARAMETERS
    ------------------------
    infile      [str]   : Name of the file containing the observation log, as written by FindHST(). 
    fields      [list]  : (Optional) If RGB fields are already known, can read in and read directly
                          in the format [R, G, B] (or [I, V, B]).
    source_id   [str]   : (Optional) Given ID of the source of interest.
    
    RETURNS
    ------------------------
    Plots a color image and returns the plt and ax objects associated with it (for further manipulation)

	"""
	
	# Reading in the observations file
    Obslog = pd.read_csv(infile).drop("Unnamed: 0", axis=1)
    
    # Finding only observations of the desired source, used to narrow down 
    # the list of observations, optional
    if source_id: 
        Obslog = Obslog[Obslog["SourceID"] == source_id].reset_index().drop("index", axis=1)
    
    if not fields: 
        display(Obslog)
        fields = input("Indices of correct RGB filters (comma, no space). If no filter available, input as 'None':").split(",")
        fields = [int(value) if value != "None" else value for value in fields]
        
        

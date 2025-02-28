(chap:astrometry)=
# Astrometric Corrections

Because we are combining the data from two different instruments, there are corrections that must be made to the coordinates of the X-ray sources such that it is propely in line with the coordinate system of *HST*. There will be some intrinsic amount of uncertainty in the XRBs' coordinates due to the differences in the X-ray and optical resolution of their respective instruments. These are accounted for by defining 1- and 2-$\sigma$ regions around each XRB, which trace out the regions within which we are 68% and 95% sure that the source falls. These regions are calculated using the information obtained through the astrometric correction process --- namely, the standard deviation of the median offset between optical and X-ray coordinates. This Chapter outlines how to conduct an astrometric correction and calculate the 1- and 2-$\sigma$ positional uncertainty regions around each XRB. 

## Calibrator Selection
To calibrate between X-ray and optical source coordinates, we need to select a sample of optical sources that we are reasonably certain are producing the X-ray emissions we detect. These include background galaxies (AGN and quasars), foreground stars, and isolated globular clusters (for example, {numref}`fig-astrometric-calibrators`). To find these sources, we must plot the coordinates of each X-ray source onto the *HST* image. We then identify the X-ray sources and their optical counterparts that are best suited for the assessment. You will want to pull the coordinates for each selected X-ray source, and (if you choose to do a by-hand correction) their associated optical source. 

```{figure} astrometric_calibrators.png
:name: fig-astrometric-calibrators
:width: 750px

Examples of X-ray sources that can be used as calibrators for the astrometric correction. The first 3 images are of AGN and quasars, based on their shapes and/or extremely red colors. The 4th is a foreground star, identifiable by the diffraction spikes. The last is a (somewhat) isolated globular cluster, based on the shape, color, and size relative to typical stars in M101.
```

## Correction Calculation and Application
As you may have guessed, I've written a function that will calculate the astrometric corrections automatically, called `CorrectAstrometry` (Chapter \ref{sec:script-correctastrometry}). It takes in the coordinates of your selected X-ray sources (`cat_coords`) and finds the nearest optical source from the coordinates of some base catalog (`base_coords`, which in this case are the coordinates from the `DAOStarFinder` catalog[^1]). Alternatively, you can calculate the astrometric correction yourself by hand-selecting the best optical source position of each calibrator and finding the median x- and y-coordinate offset between the optical and X-ray positions, as well as the standard deviation of those offsets. 

If the standard deviations on the offsets are larger than the median offsets, then that suggests the *CXO* sources are already well-aligned with the *HST* image and no coordinate shift is necessary. The standard deviations on those shifts, however, are still needed to calculate the positional uncertainties, as it represents the variations in where the optical counterparts fall with respect to the X-ray coordinates of the selected sources. Ideally these will be small, but they're expected to be non-zero. 

Below is the code I use to find the X-ray calibrators and call `CorrectAstrometry` to calculate the correction: 

```
from Sources import LoadSources
%run Sources.py
from WriteScript import WriteReg
%run WriteScript.py
from DataFrameMod import RemoveElse, FindUnique
%run DataFrameMod.py
from Align import CorrectAstrometry
%run Align.py

# Reading in sources from the DataFrame containing all of the X-ray sources
# You can use pd.read_csv() instead of LoadSources, if you prefer
M101_CSC = LoadSources("M101_CSC.frame")

# Find X-ray sources with a unique CSC ID,
# in case multiple rows in the the DataFrame have the same CSC ID
M101_unique = FindUnique(M101_CSC, header="CSC ID")

# Saving M101 to a region file that can be opened in DS9
WriteReg(M101_unique, outfile="M101_cscsources.reg", idname="CSC ID", radius=50, width=2, 
         color="hotpink", showlabel=True)

# List of source names to use as calibrators in the Astrometric Correction
# These I identified manually by inspecting the HST image in DS9 with the 
# regions saved above plotted over the image
M101_calibrators = ['2CXO J140251.0+542420', 
                    '2CXO J140248.0+542403', 
                    '2CXO J140356.0+542057',
                    '2CXO J140355.8+542058',
                    '2CXO J140357.6+541856',
                    '2CXO J140339.3+541827',
                    '2CXO J140346.1+541615',
                    '2CXO J140253.3+541855',
                    '2CXO J140252.1+541946']

print(len(M101_calibrators), "calibrators to match...")

# Using DataFrameMod.RemoveElse() to remove all but the sources above from the DataFrame
M101_calibrators = RemoveElse(M101_unique, keep=M101_calibrators, header="CSC ID")
print(len(M101_calibrators), "calibrators found.")

# Saving these as a region file, in case we want to double-check
WriteReg(M101_calibrators, outfile="M101_calibrators", radius=25, width=2)

# Setting up the base and the catalog coordinates for CorrectAstrometry
base_coords = GetCoords("M101_daofind_f555w_acs_fk5.reg")
cat_coords = [M101_calibrators["RA"].values.tolist(), M101_calibrators["Dec"].values.tolist()]

# Running CorrectAstrometry
CorrectAstrometry(base_coords, cat_coords, returnshifts=True, savebasereg="M101_astrocorrect.reg")
```


[^1]: If you had to apply a coordinate shift to the `DAOStarFinder` region file to get the `DAOStarFinder` sources to appear properly aligned to the image, then you'll want to use these shifted coordinates as the base coordinates. That makes it easily to compare the X-ray coordinates to the optical, since they'll both be visually aligned to the same image coordinate system.
# Overview of Process

For my work, my main goal was to identify XRBs within spiral galaxies, categorize them by the masses of their donor stars as low-mass, intermediate-mass, and high-mass XRBs (LMXBs, IMXBs, and HMXBs respectively), and analyze their abundances, environments, etc. The procedures I use to accomplish this are theoretically simple, but may be technically complex depending on the nuances of the particular galaxy of interest. Below I've provided a [template](https://colab.research.google.com/drive/1EXjGnlxIeJWav4ULa2UpFkykWi-bc4JP?usp=sharing)  `.ipynb` of my process (NOTE: this is currently out of date, but I'll update it soon) that works out of Google Collab (though it requires a few additional custom scripts, which are found in the `pyfiles` [of my GitHub XRBID repository](https://github.com/QosmicQi/XRBID/tree/main/xrbid_guide)). That code will not run straight out of the box on your computer without a lot of modification, but can be used to get an idea of what I'll be doing for the rest of this guide. 
In general, the process I use can be broken down into these basic steps: 

1. Select a galaxy (or galaxies) of interest;
2. Obtain the *CXO* X-ray Data from the Chandra Source Catalog (`CSC`);
3. If necessary, obtain the optical imaging from the *HST* from the Hubble Legacy Archive (HLA) or MAST either directly or through an HLA query;
4. If necessary, combine the *HST* images into a single mosaic using `AstroDrizzle` (this step isn't strictly required, but may be useful for a variety of reasons);
5. Identify point sources in HST image(s);
6. Calculate the photometric corrections on each *HST* field/image;
7. Perform astrometric corrections on the X-ray coordinates to align them to the optical image;
8. Identify candidate optical counterparts for each X-ray source;
9. Identify supernova remnants, background galaxies, and foreground stars among the optical counterparts; 
10. Extract the photometry of candidate donor stars and host star clusters, including photometric corrections; 
11. Create a color-magnitude diagram of candidate donor stars or color-color diagrams of candidate host star clusters to estimate the masses of the XRB donor stars; and
12. Perform any other analysis necessary on the XRB populations. 

## Helpful Links and Other Resources
* [Chandra Source Catalog Quick Search](http://cda.cfa.harvard.edu/cscweb/index.do)
* [The Hubble Legacy Archive](https://hla.stsci.edu/hlaview.html#)
* [The MAST database](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html) (for finding all space telescope observations of a given object), or [for Hubble only](https://mast.stsci.edu/search/ui/#/hst)
* The *HST* [WFC3 Data Handbook](https://hst-docs.stsci.edu/wfc3dhb) (for understanding how to use the data)
* A guide for *HST* [ACS data analysis](https://www.stsci.edu/hst/instrumentation/acs/data-analysis)
* [A useful guide](https://hst-docs.stsci.edu/hstdhb/4-hst-data-analysis/4-6-analyzing-hst-images) for analyzing *HST* data (general)
* *HST* [Notebook Repository](https://spacetelescope.github.io/hst_notebooks/)
* [NEARGALCAT](https://heasarc.gsfc.nasa.gov/W3Browse/galaxy-catalog/neargalcat.html) (for finding nearby galaxies to analyze; click `Browse` to use the web interface)
* The [STScI Digitized Sky Survey](https://stdatu.stsci.edu/cgi-bin/dss_form) (for pulling optical images of the galaxy)
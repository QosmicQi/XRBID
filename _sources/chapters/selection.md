(chap:selection)=
# Selecting Optical Counterparts

From this point on, the methodolgy for identifying the best candidate counterpart for each XRB and becomes fairly open ended. While there are still specific steps that need to be taken, there are a countless number of ways to do so! Here, I detail the steps I typically take to optimize my workflow. The details of ones workflow do not matter as long as you accomplish the following: 

1. Identify optical counterparts that fall within 2-$\sigma$ of each X-ray source; 

2. Within this population, identify contaminants (i.e. foreground stars, background galaxies), clusters, and point-sources (i.e. candidate donor stars);

3. Extract the photometry from clusters and point-sources; 

4. Estimate the masses of point-sources and ages of clusters to determine the most likely mass classification (low, intermediate, or high) of each XRB; and

5. Flag possible supernova remnants (SNR). 

## Incorporating Other Catalogs 

One of the things you will want to keep in mind is that some of the work you may find yourself doing here may have already been done by other studies. In particular, keep an eye out for cluster and SNR catalogs that have been published by other groups. If those exist, you will want to download those catalogs and extract the coordinates of clusters/SNRs that they've identified. You should then create a DS9 region file of those sources with `WriteScript.WriteReg()`. This will become important as you attempt to classify the optical counterparts of each X-ray source. 


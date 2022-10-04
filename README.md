# Coordinate_matching


### Description

This code is created to facilitate matching different catalogs on coordinates. It is mostly a fancy wrapper around Astropy's `match_to_catalog_sky()` method. However, this code solves a few more problems and provides additional functionality to help the matching process.
In short: this codes helps to find the optimal threshold for matching catalogs and will remove matches that are separate more than this value.
It also has the option to remove duplicates from our matching process and the process is done iteratively, to make sure no matches are lost. 

An detailed example can be found in the Jupyter Notebook. It shows how to match MaNGA to Galaxy Zoo DECaLS, based on coordinates. 
I obtained the MaNGA coordinates from [MaNGA](https://www.sdss.org/dr17/manga/) drpall, DR17. 
The Galaxy Zoo DECaLS coordinates can be found [here](https://arxiv.org/abs/2102.08414). Both were taken on 4 Oct 2022.

This code has been tested with python 3.7.9, astropy 4.2, pandas 1.2.0, numpy 1.19.2, matplotlib 3.3.2.


### Calling sequence of main functions

There are two main functions. The first helps to find the ideal separation limit that best suits the two catalogs:

```plot_separation_limit(catalog_A, catalog_B, xmin = 0.1, xmax = -1, n = 100, logscale=False)```

INPUTS  
catalog_A (list): Must be a list of two lists. First list is the ra of the first catalog, second is dec. Units for ra/dec must be deg, not sexagesimal!  
catalog_B (list): Same as catalog_A, but for the second catalog.  

OPTIONAL INPUTS  
xmin (float): Minimum separation limit to be tested, in arcsec.  
xmax (float): Maximum separation limit to be tested, in arcsec.  
n (int): Amount of points to be sampled between xmin and xmax.  
logscale (bool): Whether to plot axes in logscale.   

OUTPUTS  
None  



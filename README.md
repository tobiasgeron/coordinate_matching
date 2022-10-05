# Coordinate_matching


### Description

This code is created to facilitate matching different catalogs on coordinates. It is mostly a fancy wrapper around Astropy's `match_to_catalog_sky()` method. However, this code solves a few more problems and provides additional functionality to help the matching process.
In short: this codes helps to find the optimal separation limit for matching catalogs and will remove matches that are separated more than this value.
It also has the option to remove duplicates from our matching process and the process is done iteratively, to make sure no matches are lost. 

An detailed example can be found in the Jupyter Notebook. It shows how to match MaNGA to Galaxy Zoo DECaLS, based on coordinates. 
I obtained the MaNGA coordinates from [MaNGA](https://www.sdss.org/dr17/manga/) drpall, DR17. 
The Galaxy Zoo DECaLS coordinates can be found [here](https://arxiv.org/abs/2102.08414). Both were taken on 4 Oct 2022.

This code has been tested with python 3.7.9, astropy 4.2, pandas 1.2.0, numpy 1.19.2, matplotlib 3.3.2.


### Calling sequence

##### plot_separation_limit

There are two main functions. The first helps to find the ideal separation limit that best suits the two catalogs:

```plot_separation_limit(catalog_A, catalog_B, xmin = 0.1, xmax = -1, n = 100, logscale=False)```

INPUTS  
catalog_A (list): Must be a list of two lists. First list is the ra of the first catalog, second is dec. Units for ra/dec must be deg, not sexagesimal!  
catalog_B (list): Same as catalog_A, but for the second catalog.  

OPTIONAL INPUTS  
xmin (float): Minimum separation limit to be tested, in arcsec. Default is 0.1 arcsec.  
xmax (float): Maximum separation limit to be tested, in arcsec. If left at -1, the code will choose an appropriate value automatically.  
n (int): Amount of points to be sampled between xmin and xmax.  Default is 100.  
logscale (bool): Whether to plot axes in logscale. Default is False.  

OUTPUTS  
None  

NOTES  
For speed purposes, we do not remove duplicates or perform the matching iteratively in this function. 
This is done in the more detailed `match_catalogs()` function. Therefore, the amount of matches shown is 
only an estimate and should only serve as a rough guideline. 

The output plot can be used to determine the ideal separation limit. Please refer to the example Jupyter Notebook for more details.

##### match_catalogs

After an ideal separation limit is determined, this next function actually does the matching:

```match_catalogs(catalog_A, catalog_B, remove_duplicates = True, limit = 0, recursive = True, max_loops = 10, i_loop = 0):```

DESCRIPTION  
Will map catalog_B to catalog_A (i.e. for every element in A, find closest element in B).  

INPUTS  
catalog_A (list): Must be a list of two lists. First list is the ra of the first catalog, second is dec. Units for ra/dec must be deg, not sexagesimal!  
catalog_B (list): Same as catalog_A, but for the second catalog.  

OPTIONAL INPUTS  
limit (float): The separation limit, in arcsec. All targets that are further away than this limit will be removed and replaced by np.nan. Default is 10 arcsec.  
remove_duplicates (bool): Whether to remove duplicate matches. Default and recommended setting is 'True'.  
recursive (bool): Whether to match recursively. Default and recommended setting is 'True'.  
max_loops (int): How many iterations we can do for the recursive matching. Default is 10. 
i_loop (int): Which iteration we are currently on.   

OUTPUT  
idx (list): List of length equal to catalog_A. Every entry has the index of corresponding row in catalog_B. If no match is found, entry is np.nan instead.   
d2d (list): List of length equal to catalog_A. Contains on-sky distance to match in catalog_B.   
n_removed (int): Used for recursive matching. It notes how many duplicates are removed in this iteration.  

NOTES  
The code is significantly faster if you choose the smaller catalog as catalog A!  



More details can be found in the Jupyter Notebook, as well as the `.py` file containing the code.




### TODO:
Explain what recursive matching actually does..



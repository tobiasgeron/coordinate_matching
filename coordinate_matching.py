'''
Made by Tobias Géron on 4 Oct 2022, based on older code on my laptop. 

TODO:
Make actual matching function.
Add redshift in matching process?
Automate the finding of separation limit? Can do gaussian smoothing, apply differential and see where is maximum. That'll probably be a good first approximation.
'''

###############
### IMPORTS ###
###############

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u


############
### CODE ###
############

def plot_separation_limit(catalog_A, catalog_B, xmin = 0.1, xmax = -1, n = 100, logscale=False):

    '''
    DESCRIPTION
    This function will plot a figure where you can see which separation limit is appropriate to use while matching two catalogs.
    Will match catalog_B to catalog_A (i.e. for every element in A, find closest element in B).

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
    '''


    assert len(catalog_A) == 2 and len(catalog_B) == 2, 'Catalog_A must be a list of 2 arrays. First array is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.'
    
    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]


    #set xmax depending on logscale
    if logscale:
        if xmax == -1:
            xmax = 60*60*10 #so 10 degrees
    else:
        if xmax == -1:
            xmax = 50


    #space sep_limits correctly
    if logscale:
        sep_limits = np.logspace(np.log10(xmin), np.log10(xmax), n)
    else:
        sep_limits = np.linspace(xmin, xmax, n)

        
    Catalog_A = SkyCoord(ra=catalog_A[0]*u.degree, dec=catalog_A[1]*u.degree)
    Catalog_B = SkyCoord(ra=catalog_B[0]*u.degree, dec=catalog_B[1]*u.degree)

    idx, d2d, _ = Catalog_A.match_to_catalog_sky(Catalog_B)

    n_matches = []
    for sep in sep_limits:
        n_matches.append(len(np.where(d2d < sep*u.arcsec)[0]))


    fig, ax = plt.subplots(figsize=(8.5,5))
    ax.plot(sep_limits,n_matches)

    if logscale:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.axvline(x=1,lw=1, ls = ':', c='black')
        ax.axvline(x=60,lw=1, ls = ':', c='black')
        ax.axvline(x=60*60,lw=1, ls = ':', c='black')


    ax.set_xlabel('Separation limit [arcsec]', fontsize = 16)
    ax.set_ylabel('Amount of matches found', fontsize = 16)
    plt.show()


catalog_match_plot_separation_radius = plot_separation_limit # for backwards compatibility


def match_catalogs(catalog_A, catalog_B, limit = 10., remove_duplicates = True, recursive = True, max_loops = 10, i_loop = 0):
    '''
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
    '''

    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]

    
    Catalog_A = SkyCoord(ra=catalog_A[0]*u.degree, dec=catalog_A[1]*u.degree)
    Catalog_B = SkyCoord(ra=catalog_B[0]*u.degree, dec=catalog_B[1]*u.degree)

    idx, d2d, _ = Catalog_A.match_to_catalog_sky(Catalog_B)
    idx_final = list(idx)
    d2d_arcsec = d2d.arcsec
    
    # Removing targets that are too far away,
    for i in range(len(d2d)):
        # if it is within the limit
        if d2d_arcsec[i] >= limit:
            idx_final[i] = np.nan
            d2d[i] = np.nan
    
    n_removed = 0
    #Removing duplicates, if wanted
    if remove_duplicates:
        df = pd.DataFrame(idx_final)
        idxs_dupl = np.where(df.duplicated(keep=False))[0]
        #can drop NaN here? Don't think so, cause then indexing is wrong
        
        for idx_dupl in idxs_dupl:
            if np.isnan(idx_final[idx_dupl]):
                continue
                
            idxs = np.where(np.array(idx_final) == idx_final[idx_dupl])[0] #find the similar idxs
            closest_idx = idxs[np.where(d2d[idxs] == np.min(d2d[idxs]))[0][0]] #find which ones of them is closest
            
            for j in idxs:
                if j!= closest_idx:
                    idx_final[j] = np.nan
                    d2d[j] = np.nan
                    n_removed+=1

    idx = np.array(idx_final)
    if recursive and n_removed > 0 and i_loop < max_loops:
        #Find targets in Catalog A and B without a match (maybe only do ones in n_removed?)
        idxs_A_nomatch = np.where(np.isnan(idx))[0]
        idxs_B_nomatch = [i for i in range(len(catalog_B[0])) if i not in idx] #this takes a long time if catalog B is long, speed this up?
        
        idx_temp, d2d_temp, n_removed = match_catalogs([catalog_A[0][idxs_A_nomatch], catalog_A[1][idxs_A_nomatch]],
                                                                [catalog_B[0][idxs_B_nomatch], catalog_B[1][idxs_B_nomatch]],
                                                                remove_duplicates = remove_duplicates, limit = limit, 
                                                                recursive = recursive, i_loop = i_loop+1, max_loops = max_loops)

        # Combine idx and d2d
        for i in range(len(idx)):
            if i in idxs_A_nomatch: #if it didn't have a match before
                j = np.where(idxs_A_nomatch == i)[0][0]
                if ~np.isnan(idx_temp[j]): #if it has a match in the new search
                    idx[i] = idxs_B_nomatch[int(idx_temp[j])]
                    d2d[i] = d2d_temp[j]



        
    return idx, d2d, n_removed


def merge_catalogs(catalog_A, catalog_B, idx, suffixes = ('_A','_B')):
    '''
    DESCRIPTION
    After doing match_catalogs(catalog_A, catalog_B), this function takes the output of that and actually merges the catalogs. 
    Assumes both catalogs are pandas dfs.

    INPUTS:
    catalog_A (pandas DataFrame): first catalog
    catalog_B (pandas DataFrame): second catalog
    idx: first output of match_catalogs(catalog_A, catalog_B). For every row in catalog_A, maps to a row in catalog B. This is what 
    we use to merge.

    OPTIONAL INPUTS:
    suffixes (list): if there are duplicate columns, what suffixes to use for each catalog.

    NOTES
    I'm assuming catalog_A is the smaller one (as advised above). 
    This code also assumes that there is no "index" column in either df. If so, it will be overwritten.
    
    '''

    assert len(idx) == len(catalog_A),'The length of idx should be the same as catalog A. For every row in catalog_A, idx should map to a row in catalog B. This is what we use to merge.'


    # Add index to catalog_B
    catalog_B['index'] = list(range(len(catalog_B)))

    # Remove all in catalog A where idx is np.nan, so no match was found
    ind_selection = np.where(np.isnan(idx))[0]
    catalog_A = catalog_A.drop(ind_selection,axis=0)
    catalog_A = catalog_A.reset_index(drop=True)

    # Get rid of all the NaN entries
    idx_denan = [int(i) for i in idx if ~np.isnan(i)]
    
    # Merge
    catalog_A['index'] = list(catalog_B['index'][idx_denan]) #Can just do idx_denan?
    catalog_C = catalog_A.merge(catalog_B,left_on='index',right_on='index', how='inner', suffixes = suffixes)
    catalog_C = catalog_C.drop(columns=['index'])

    return catalog_C



def plot_on_sky(catalog_A, catalog_B, labels = [1,2], markersize = 1.0, alpha = 0.3, frac_plot = 1):
    '''
    DESCRIPTION
    Takes the coordinates of both catalogs and plots them on the sky. This helps visualise whether they actually overlap.

    INPUTS
    catalog_A (list): Must be a list of two lists. First list is the ra of the first catalog, second is dec. Units for ra/dec must be deg, not sexagesimal!
    catalog_B (list): Same as catalog_A, but for the second catalog.

    OPTIONAL INPUTS
    labels (list): List of strings. Will be the labels used in subscript for the two catalogs. Default is [1,2].
    marksersize (float): Size of the markers. Default is 1.0.
    alpha (float): opacity of the markers. Default is 0.3.
    frac_plot (float): between 0 and 1. Determines fraction of catalog to actually plot. For large catalogs, 
    plotting >1M dots will take too much time, and is not useful. If you just want an idea of where the footprint is, you can reduce this.

    OUTPUT
    None

    NOTES
    I'm assuming catalog_A is the smaller one (as advised above). Therefore, I'm plotting catalog_B first, so more of catalog_A is visible.
    
    TODO
    Can add give option to only plot a fraction of the catalogs. Will speed up larger catalogs
    '''

    assert len(catalog_A) == 2 and len(catalog_B) == 2, 'Catalog_A must be a list of 2 arrays. First array is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.'
    
    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]


    # Frac plot
    if frac_plot != 1:
        inds_A = np.random.choice(range(len(catalog_A[0])), replace = False, size = int(frac_plot * len(catalog_A[0])))
        catalog_A = np.array(catalog_A).T[inds_A].T

        inds_B = np.random.choice(range(len(catalog_B[0])), replace = False, size = int(frac_plot * len(catalog_B[0])))
        catalog_B = np.array(catalog_B).T[inds_B].T
    
    
    if len(catalog_B[0]) > 0:
        plt.scatter(catalog_B[0], catalog_B[1], s=markersize, alpha=alpha, label = labels[1]) #plotting B first, as we assume B is the smalelr catalog
    if len(catalog_A[0]) > 0:
        plt.scatter(catalog_A[0], catalog_A[1], s=markersize, alpha=alpha, label = labels[0])

    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
    
    plt.legend()
    
    plt.show()


def plot_coordinate_difference(catalog_A, catalog_B, labels = [1,2]):
    '''
    DESCRIPTION
    After matching, this function helps to check how accurate the matches are.

    INPUTS
    catalog_A (list): Must be a list of two lists. First list is the ra of the first catalog, second is dec. Units for ra/dec must be deg, not sexagesimal!
    catalog_B (list): Same as catalog_A, but for the second catalog. It is implied that this is after matching the catalogs, so that the first element in catalog_B is matched to the first element in catalog_A, and so on.

    OPTIONAL INPUTS
    labels (list): List of strings. Will be the labels used in subscript for the two catalogs. Default is [1,2].

    OUTPUT
    None

    NOTES
    It is implied that this is after matching the catalogs, so that the first element in catalog_B is matched to the first element in catalog_A, and so on.
    
    '''

    assert len(catalog_A) == 2 and len(catalog_B) == 2, 'Catalog_A must be a list of 2 arrays. First array is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.'
    
    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]

    # Can just use astropy functions here to be more acurate? Shouldn't make much of a difference.
    delta_ra = (catalog_A[0] - catalog_B[0]) * 3600
    delta_dec = (catalog_A[1] - catalog_B[1]) * 3600
    delta = np.sqrt(delta_ra**2 + delta_dec**2)

    
    plt.figure(figsize = (10,5))
    plt.subplot(1,2,1)
    plt.scatter(delta_ra,delta_dec)
    plt.xlabel(r'RA$_{' + str(labels[0]) + r'}$ - RA$_{' + str(labels[1]) + '}$ [arcsec]', fontsize = 14)
    plt.ylabel(r'DEC$_{' + str(labels[0]) + r'}$ - DEC$_{' + str(labels[1]) + '}$ [arcsec]', fontsize = 14)
    

    plt.subplot(1,2,2)
    plt.hist(delta)
    plt.xlabel(r'$\sqrt{ \Delta RA + \Delta DEC}$  [arcsec]', fontsize = 14)
    plt.ylabel('Frequency', fontsize = 14)

    plt.tight_layout()
    plt.show()


def sexagesimal_to_degree(ra,dec):
    '''
    DESCRIPTION
    Using astropy to transform coordinates from sexagesimal (in hms for ra and dms for dec) to degrees.

    INPUTS
    ra (str): ra in sexagesimal notation
    dec (str): dec in sexagesimal notation

    OPTIONAL INPUTS
    

    OUTPUT
    ra (float): ra in degrees
    dec (float): dec in degrees

    NOTES
    
    '''
    c = SkyCoord(ra = ra, dec = dec,unit=(u.hourangle, u.deg))
    
    return c.ra.degree, c.dec.degree
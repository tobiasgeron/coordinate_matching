'''
Made by Tobias Géron on 4 Oct 2022, based on older code on my laptop. 

TODO:
Make plotting function for delta ra/dec
Make actual matchting function.
'''

###############
### IMPORTS ###
###############

import numpy as np
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
import pandas as pd


############
### CODE ###
############


def plot_separation_radius(catalog_A, catalog_B, xmin = 0.1, xmax = -1, n = 100, logscale=False):

    '''
    DESCRIPTION
    This function will plot a figure where you can see which separation limit is appropriate to use while matching two catalogs.
    Will match catalog_B to catalog_A (i.e. for every element in A, find closest element in B).

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

    NOTES
    For speed purposes, we do not removed duplicates or perform the matching iteratively in this function. 
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

    idx, d2d, d3d = Catalog_A.match_to_catalog_sky(Catalog_B)

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


catalog_match_plot_separation_radius = plot_separation_radius # for backwards compatibility


def match_catalogs(catalog_A, catalog_B, remove_duplicates = True, limit = 0, recursive = True, 
                                max_loops = 10, i_loop = 0):
    '''
    If limit != 0, will remove all targets who are separated more than that limit (in arcsec)
    If remove_duplicates == True, will only keep the closest target
    Will match catalog_B to catalog_A (i.e. for every element in A, find closest element in B)
    
    Significantly faster to take the smaller catalog as catalog A
    
    Take care when debugging, a lot of different kinds of indices... Confusing.
    '''

    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]

    
    Catalog_A = SkyCoord(ra=catalog_A[0]*u.degree, dec=catalog_A[1]*u.degree)
    Catalog_B = SkyCoord(ra=catalog_B[0]*u.degree, dec=catalog_B[1]*u.degree)

    idx, d2d, d3d = Catalog_A.match_to_catalog_sky(Catalog_B)
    idx_final = list(idx)
    d2d_arcsec = d2d.arcsec
    
    # Removing targets that are too far away, if specified
    if limit > 0:
        for i in range(len(d2d)):
            # if it is within the limit
            if d2d_arcsec[i] >= limit:
                idx_final[i] = np.nan
                d2d[i] = np.nan
                d3d[i] = np.nan
    
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
                    d3d[j] = np.nan
                    n_removed+=1

    idx = np.array(idx_final)
    if recursive and n_removed > 0 and i_loop < max_loops:
        #Find targets in Catalog A and B without a match (maybe only do ones in n_removed?)
        idxs_A_nomatch = np.where(np.isnan(idx))[0]
        idxs_B_nomatch = [i for i in range(len(catalog_B[0])) if i not in idx] #this takes a long time if catalog B is long, speed this up?
        
        idx_temp, d2d_temp, d3d_temp, n_removed = match_catalogs([catalog_A[0][idxs_A_nomatch], catalog_A[1][idxs_A_nomatch]],
                                                                [catalog_B[0][idxs_B_nomatch], catalog_B[1][idxs_B_nomatch]],
                                                                remove_duplicates = remove_duplicates, limit = limit, 
                                                                recursive = recursive, i_loop = i_loop+1, max_loops = max_loops)

        # Combine idx, d2d and d3d
        for i in range(len(idx)):
            if i in idxs_A_nomatch: #if it didn't have a match before
                j = np.where(idxs_A_nomatch == i)[0][0]
                if ~np.isnan(idx_temp[j]): #if it has a match in the new search
                    idx[i] = idxs_B_nomatch[int(idx_temp[j])]
                    d2d[i] = d2d_temp[j]
                    d3d[i] = d3d_temp[j]



        
    return idx, d2d, d3d, n_removed


def plot_on_sky(catalog_A, catalog_B, labels = [1,2], markersize = 1, alpha = 0.3):
    assert len(catalog_A) == 2 and len(catalog_B) == 2, 'Catalog_A must be a list of 2 arrays. First array is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.'
    
    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]

    plt.scatter(catalog_B[0], catalog_B[1], s=markersize, alpha=alpha, label = labels[1]) #plotting B first, as we assume B is the smalelr catalog
    plt.scatter(catalog_A[0], catalog_A[1], s=markersize, alpha=alpha, label = labels[0])

    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
    
    plt.legend()
    
    plt.show()



def plot_coordinate_difference(catalog_A, catalog_B, labels = [1,2]):
    '''
    '''

    assert len(catalog_A) == 2 and len(catalog_B) == 2, 'Catalog_A must be a list of 2 arrays. First array is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.'
    
    # Transform to np.array
    catalog_A = [np.array(catalog_A[0]), np.array(catalog_A[1])]
    catalog_B = [np.array(catalog_B[0]), np.array(catalog_B[1])]

    delta_ra = (catalog_A[0] - catalog_B[0]) * 3600
    delta_dec = (catalog_A[1] - catalog_B[1]) * 3600
    delta = np.sqrt(delta_ra**2 + delta_dec**2)

    
    plt.figure(figsize = (10,5))
    plt.subplot(1,2,1)
    plt.scatter(delta_ra,delta_dec)
    plt.xlabel(r'RA$_{ \rm' + str(labels[0]) + r'}$ - RA$_{ \rm' + str(labels[1]) + '}$ [arcsec]', fontsize = 14)
    plt.ylabel(r'DEC$_{ \rm' + str(labels[0]) + r'}$ - DEC$_{ \rm' + str(labels[1]) + '}$ [arcsec]', fontsize = 14)
    

    plt.subplot(1,2,2)
    plt.hist(delta)
    plt.xlabel(r'$\sqrt{ \Delta RA + \Delta DEC}$  [arcsec]', fontsize = 14)
    plt.ylabel('Frequency', fontsize = 14)

    plt.tight_layout()
    plt.show()







def sexagesimal_to_degree(ra,dec):
    '''
    Using astropy to transform coordinates from sexagesimal (in hms for ra and dms for dec) to degrees
    '''
    c = SkyCoord(ra = ra, dec = dec,unit=(u.hourangle, u.deg))
    
    return c.ra.degree, c.dec.degree



import numpy as np
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
import pandas as pd
from tqdm.notebook import tqdm



def catalog_match_plot_separation_radius(catalog_A, catalog_B, xmin = -1, xmax = -1, n = 100, logscale=False):

    '''
    This function will plot an easy figure where you can see which separation limit is appropriate for you two catalogs. 
    Catalog_A must be a list of 2 arrays. First element is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.
    Will match catalog_B to catalog_A (i.e. for every element in A, find closest element in B)
    
    xmin and xmax are the minimum and maximum separation limits to be tested, in arcsec.
    n is the amount of datapoints that should be tested. 

    logscale will put both axes on logscale. Recommended only if you want to probe until degree range. If turned on, will probe between 1 arcsec and 10 degs.
    '''


    assert len(catalog_A) == 2 and len(catalog_B) == 2, 'Catalog_A must be a list of 2 arrays. First element is ra, second is dec. Idem for catalog_B. Units for ra/dec must be deg.'

    #set limits of sep_limit correctly
    if logscale:
        if xmax == -1:
            xmax = 60*60*10 #so 10 degrees
    else:
        if xmax == -1:
            xmax = 50
    if xmin == -1:
        xmin = 0.1


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


def match_catalogs(catalog_A, catalog_B, remove_duplicates = True, limit = 0, recursive = True, 
                                max_loops = 10, i_loop = 0):
    '''
    If limit != 0, will remove all targets who are separated more than that limit (in arcsec)
    If remove_duplicates == True, will only keep the closest target
    Will match catalog_B to catalog_A (i.e. for every element in A, find closest element in B)
    
    Faster to take the smaller catalog as catalog A
    
    Take care when debugging, a lot of different kinds of indices... Confusing.
    '''
    
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
        idxs_B_nomatch = [i for i in tqdm(range(len(catalog_B[0]))) if i not in idx] #this takes a long time if catalog B is long, speed this up?
        
        idx_temp, d2d_temp, d3d_temp, n_removed = match_catalogs_recursively([catalog_A[0][idxs_A_nomatch], catalog_A[1][idxs_A_nomatch]],
                                                                [catalog_B[0][idxs_B_nomatch], catalog_B[1][idxs_B_nomatch]],
                                                                remove_duplicates = remove_duplicates, limit = limit, 
                                                                recursive = recursive, i_loop = i_loop+1, max_loops = max_loops)

        # Combine idx, d2d and d3d
        for i in tqdm(range(len(idx))):
            if i in idxs_A_nomatch: #if it didn't have a match before
                j = np.where(idxs_A_nomatch == i)[0][0]
                if ~np.isnan(idx_temp[j]): #if it has a match in the new search
                    idx[i] = idxs_B_nomatch[int(idx_temp[j])]
                    d2d[i] = d2d_temp[j]
                    d3d[i] = d3d_temp[j]



        
    return idx, d2d, d3d, n_removed
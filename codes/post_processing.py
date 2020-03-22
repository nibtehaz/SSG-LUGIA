"""
    Performs the necessary post-processing
"""

import numpy as np

def medianFiltering(ys, model_params):
    """
    Applies median filter on the computed distance values
    
    Arguments:
        ys {numpy.array} -- array of distance values
        model_params {dictionary} -- SSG-LUGIA model configuration
    
    Returns:
        numpy.array -- median filtered array of distance values
    """

    
    ys2 = np.array(ys)
    
    wlen = model_params["median_filter_window_len"]
    wlen_2 = wlen // 2                          # computing half length for efficient computation

    for i in range(wlen_2,len(ys)-wlen_2):

            ys2[i] = np.median(ys[i-wlen_2:i+wlen_2])       # applying median filter
            
    return ys2


def binaizeFilteredDistance(ys_mf):
    """
    Binarizing by putting threshold on the filtered distance 
    
    Arguments:
        ys_mf {numpy array} -- filtered distance array

    Returns:
        numpy.array -- thresholded distance array
    """

    yp_mf = np.array(ys_mf)                 # initialize an array
    yp_mf[np.where(ys_mf<0)] = -1           # -1 => Native 
    yp_mf[np.where(ys_mf>=0)] = 1           # +1 => Alien

    return yp_mf


def inferLabel(y, genome_length, model_params):
    """
    Infers the labels of the individual nucleotides from the
    prediction of the windowed sequences
    
    Arguments:
        y {numpy array} -- array containing predictions for windowed sequences
        genome_length {int} -- length of the genome
        model_params {dictionary} -- SSG-LUGIA model configuration
    
    Returns:
        numpy.array -- array containing labels for individual nucleotides
    """

    step_size = model_params["dw"]          # step size in moving window computation

    Y = np.zeros((genome_length,1))         # initialization

    for i in range(len(y)):

        if(y[i]==-1):
            Y[i*step_size:(i+1)*step_size] = 1      # updating based on step size
    

    return Y[:]


def trimRegions(y, model_params):
    """
    Removes smaller islands
    
    Arguments:
        y {numpy array} -- prediction for the entire genome
        model_params {dictionary} -- SSG-LUGIA model configuration
    
    Returns:
        numpy.array -- processed prediction for the entire genome
    """


    st_island = -1              

    for i in (range(len(y))):

        pt = y[i][0]                # value at point

        if(pt>0):
            if(st_island==-1):
                st_island = i       # initialize
            else:
                pass
        
        else:
            if(st_island != -1):
                en_island = i-1        # update

                island_len = en_island - st_island + 1

                if not(island_len>=model_params['min_island_len']):
                                    
                    for j in range(st_island,en_island+1):
                        y[j][0] = 0

                st_island = -1          # re-initialize

            else:
                pass
                                        # final case
    if(st_island != -1):

        en_island = len(y)-1

        island_len = en_island - st_island + 1

        if not (island_len>=model_params['min_island_len']):
            
            for j in range(st_island,en_island+1):
                y[j][0] = 0

    
    return y
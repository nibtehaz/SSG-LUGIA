"""
    Wrapper script of the SSG-LUGIA pipeline
"""

import numpy as np
from file_handling import loadGenome
from models import loadModel, loadModelJson
from feature_extraction import extractFeatures
from anomaly_detection import detectAnomalies
from post_processing import medianFiltering, binaizeFilteredDistance, inferLabel, trimRegions
from utils import inputModel, getGenomicIslands


def SSG_LUGIA(sequence_fasta_file_path = None, genome_sequence = None, model_name = None, model_parameters = None):
    """
    SSG-LUGIA Pipeline
    
    Keyword Arguments:
        sequence_fasta_file_path {string} -- path to the genome sequence fasta file (default: {None})
        genome_sequence {string} -- the actual genome sequence (default: {None})
        model_name {string} -- name of one of the standard models or path to model json file (default: {None})
        model_parameters {dictionary} -- SSG-LUGIA model configuration (default: {None})
    
    Returns:
        list of tuples : list of the tuples containing genomic islands start and end coordinates
    """

    if(sequence_fasta_file_path != None):               # if a genome sequence fasta file has been provided

        genome = loadGenome(sequence_fasta_file_path)
    
    elif(genome_sequence != None):                      # otherwise the genome sequence has been provided

        genome = genome_sequence

    else:                                               # no genome sequence has been provided

        raise ValueError('Please input either genome sequence of sequence fasta file path')


    if(model_name != None):                             
                                                # if a standard model name has been provided
        if (model_name in ['SSG-LUGIA-F','SSG-LUGIA-R','SSG-LUGIA-P']):
            model = loadModel(model_name)
                                                # if path to the json file has been provided
        else:
            model = loadModelJson(model_name)            
    
    elif(model_parameters != None):                     # a custom model parameter configuration has been provided

        model = model_parameters

    else:                                               # interactively input the model

        model = inputModel()

    

    X = extractFeatures(genome,model)                   # extracting features

    print()

    print('Detecting Anomalies')

    (yp,ys) = detectAnomalies(X, model)                 # Detecting anomalies

    print('Post-Processing')

    ys_mf = medianFiltering(ys,model)                   # applying median filter
    
    yp_mf = binaizeFilteredDistance(ys_mf)              # binarizing the values

    label_pred = inferLabel(yp_mf, len(genome), model)  # infer prediction

    label = trimRegions(label_pred, model)              # eliminate small islands

    genomic_islands = getGenomicIslands(label)          # extract genomic islands

    print('Total Number of Genomic Islands = {}'.format(len(genomic_islands)))

    print('Genomic Islands :')
    
    for genomic_island in genomic_islands:
        print('{}\t{}'.format(genomic_island[0], genomic_island[1]))

    return genomic_islands                              # return the islands



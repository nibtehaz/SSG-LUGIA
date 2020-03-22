import numpy as np
from file_handling import loadGenome
from models import loadModel
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
        model_name {string} -- name of one of the standard models (default: {None})
        model_parameters {dictionary} -- SSG-LUGIA model configuration (default: {None})
    
    Returns:
        list of tuples : list of the tuples containing genomic islands start and end coordinates
    """

    if(sequence_fasta_file_path != None):

        genome = loadGenome(sequence_fasta_file_path)
    
    elif(genome_sequence != None):

        genome = genome_sequence

    else:

        raise ValueError('Please input either genome sequence of sequence fasta file path')


    if(model_name != None):

        model = loadModel(model_name)
    
    elif(model_parameters != None):

        model = model_parameters

    else:

        model = inputModel()

    

    X = extractFeatures(genome,model)

    print()

    print('Detecting Anomalies')

    (yp,ys) = detectAnomalies(X, model)

    print('Post-Processing')

    ys_mf = medianFiltering(ys,model)
    
    yp_mf = binaizeFilteredDistance(ys_mf)

    label_pred = inferLabel(yp_mf, len(genome), model)

    label = trimRegions(label_pred, model)

    genomic_islands = getGenomicIslands(label)

    print('Total Number of Genomic Islands = {}'.format(len(genomic_islands)))

    print('Genomic Islands :')
    
    for genomic_island in genomic_islands:
        print('{}\t{}'.format(genomic_island[0], genomic_island[1]))

    return genomic_islands


#Execute(sequence_fasta_file_path='sample_data/NC_003198.1.fasta',model_name='SSG-LUGIA-F')

"""
    Utility functions
"""

def inputModel():
    """
    Provides a command line interface to input the model

    Returns:
        dictionary -- SSG-LUGIA model configuration
    """

    print('--------------------')
    print('|Feature Extraction|')
    print('--------------------')
    

    w = int(input('>Length of sliding window = '))
    dw = int(input('>Step size of sliding window = '))

    while(True):
        karlin_mode = int(input('>Mode of Karlin Bias : [1/2/3]\n[1]-normalized\n[2]-original\n[3]-plain 2-mer frequency\n>'))

        if(karlin_mode==1):
            karlin_mode = "normalized"
            print("Karlin Bias mode selected = 'normalized'")
            break

        
        elif(karlin_mode==2):
            karlin_mode = "original"
            print("Karlin Bias mode selected = 'original'")
            break

        elif(karlin_mode==3):
            karlin_mode = "raw"
            print("Karlin Bias mode selected = 'raw'")
            break
        
        else:
            print("Please select an option from 1 or 2 or 3")

    while(True):

        pca_dn = int(input('>No. of PCA Components for Dinucleotide Bias (values : 1~16) = '))

        if(pca_dn<1 or pca_dn>16):
            print("Please select a value from 1 to 16")
        else:
            break

    while(True):

        pca_amino_acid = int(input('>No. of PCA Components for Codon Bias (values : 1~20) = '))

        if(pca_amino_acid<1 or pca_amino_acid>20):
            print("Please select a value from 1 to 20")
        else:
            break
    
    while(True):

        pca_kmer4 = int(input('>No. of PCA Components for 4-mer Frequency (values : 1~256) = '))

        if(pca_kmer4<1 or pca_kmer4>256):
            print("Please select a value from 1 to 256")
        else:
            break
    
    while(True):
        entropy_features = input('>Use Entropy based features : [y/n]\n[y]-yes\n[n]-no\n>')

        if(entropy_features=='y'):
            entropy_features = True
            print("Entropy based features included")
            break

        
        elif(entropy_features=='n'):
            entropy_features = False
            print("Entropy based features discarded")
            break

        else:
            print("Please select an option from y or n")

    
    print('-------------------')
    print('|Anomaly Detection|')
    print('-------------------')

    contamination_model1 = float(input('>Contamination for first level model (floating number from 0 to 0.99) = '))

    support_fraction_model1 = 0.75

    contamination_model2 = float(input('>Contamination for second level model (floating number from 0 to 0.99) = '))

    support_fraction_model2 = 0.9

    print('----------------')
    print('|Postprocessing|')
    print('----------------')
    median_filter_window_len = int(input('>Length of median filter window = '))
    min_island_len = int(input('>Minimum length of genomic island = '))

    model_params = {    
                "w" : w,
                "dw" : dw,                
                "karlin_mode" : karlin_mode,
                "pca_dn" : pca_dn,
                "pca_amino_acid" : pca_amino_acid,
                "pca_kmer4" : pca_kmer4,
                "entropy_features" : entropy_features,
                "contamination_model1" : contamination_model1,
                "support_fraction_model1" : support_fraction_model1,
                "contamination_model2" : contamination_model2,
                "support_fraction_model2" : support_fraction_model2,
                "median_filter_window_len" : median_filter_window_len,
                "min_island_len" : min_island_len
            }

    return model_params


def getGenomicIslands(label):
    """
    Extracts the genomics islands from the label
    
    Arguments:
        label {list} -- predicted label

    Returns:
        list of tuples -- list of genomic islands
    """    

    gi_islands = []                 # genomic islands

    st_island = -1
    

    for i in (range(len(label))):

        pt = label[i][0]

        if(pt>0):
            if(st_island==-1):
                st_island = i
            else:
                pass
        
        else:
            if(st_island != -1):
                en_island = i-1                
            
                gi_islands.append((st_island+1,en_island))
                                    
                st_island = -1

            else:
                pass

    if(st_island != -1):
        en_island = len(label)-1
    
        gi_islands.append((st_island+1,en_island))
            
   
    return gi_islands
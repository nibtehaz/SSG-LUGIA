import json

"""
    Trained models presented in the manuscript.
    Provides method to load the models
"""

SSG_LUGIA_F = {                
                "w":10000,
                "dw":100,                
                "karlin_mode":"normalized",
                "pca_dn":2,
                "pca_amino_acid":2,
                "pca_kmer4":2,
                "entropy_features":True,                
                "contamination_model1":0.15,
                "support_fraction_model1":0.75,
                "contamination_model2":0.05,                
                "support_fraction_model2":0.9,
                "median_filter_window_len":400,
                "min_island_len":10000
            }


SSG_LUGIA_R = {    
                "w":10000,
                "dw":100,                
                "karlin_mode":"normalized",
                "pca_dn":2,
                "pca_amino_acid":2,
                "pca_kmer4":2,
                "entropy_features":True,                
                "contamination_model1":0.2,
                "support_fraction_model1":0.75,
                "contamination_model2":0.25,
                "support_fraction_model2":0.9,
                "median_filter_window_len":400,
                "min_island_len":10000
            }

SSG_LUGIA_P = {    
                "w":10000,
                "dw":100,
                "karlin_mode":"normalized",
                "pca_dn":2,
                "pca_amino_acid":2,
                "pca_kmer4":2,
                "entropy_features":True,                
                "contamination_model1":0.075,
                "support_fraction_model1":0.75,
                "contamination_model2":0.075,                
                "support_fraction_model2":0.9,
                "median_filter_window_len":400,                
                "min_island_len":10000
          }

def loadModel(model_name):
    """
    Loads one of the standard models
    
    Arguments:
        model_name {string} -- name of the model

    Returns:
        dictionary -- the pre-tuned model
    """

    if(model_name=='SSG-LUGIA-F'):
        return SSG_LUGIA_F
    
    elif(model_name=='SSG-LUGIA-R'):
        return SSG_LUGIA_R

    elif(model_name=='SSG-LUGIA-P'):
        return SSG_LUGIA_P

    else:
        raise ValueError('Invalid Model Name Provided')


def loadModelJson(model_path):
    """
    Loads a SSG-LUGIA model created by the user saved as .json file
    
    Arguments:
        model_path {string} -- path to the model json file

    Returns:
        dictionary -- the user-defined model
    """

    try:
        fp = open(model_path,'r')
        model_data = fp.read()
        fp.close()

        try:
            model = json.loads()

            return model

        except:
            raise ValueError('Invalid Json File Provided')
    
    except:
        raise ValueError('Invalid Path Provided')
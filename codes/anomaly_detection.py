"""
    Detects the anomalies based on the computed features
"""

from sklearn.covariance import EllipticEnvelope
from sklearn.metrics import classification_report
import numpy as np
import warnings
warnings.filterwarnings("ignore")           # ignores the warnings raised by scikit-learn


def detectAnomalies(X,model_params):
    """
    Detects the anomalies using Mahalonobis Distance
    
    Arguments:
        X {2d numpy.array} -- features of the windowed sequences
        model_params {dictionary} -- SSG-LUGIA model configuration    
    
    Returns:
        yp {numpy.array} -- binary prediction of the anomalies
        ys {numpy.array} -- mahalonobis distance of the anomalies
    """


                # we use the EllipticEnvelope model from Scikit-Learn library
                # to detect anomalies using mahalonobis distance    
    elenv = EllipticEnvelope(contamination=model_params['contamination_model1'],support_fraction=model_params['support_fraction_model1'],random_state=3)

    elenv.fit(X)

    yp = elenv.predict(X)                   # binary prediction
    ys = elenv.decision_function(X)         # mahalonobis distance computation

    

    X2 = X[np.where(yp==1)]                 # selecting only the windows predicted native


                                            # performing anomaly detection again
    elenv2 = EllipticEnvelope(contamination=model_params['contamination_model2'],support_fraction=model_params['support_fraction_model2'],random_state=3)
    elenv2.fit(X2)

    yp2 = elenv2.predict(X2)                   # binary prediction
    ys2 = elenv2.decision_function(X2)         # mahalonobis distance computation

            
    ys[np.where(yp==1)] = ys2               # updating the binary prediction based on level 2 detection
    yp[np.where(yp==1)] = yp2               # updating the mahalonobis distance based on level 2 detection


    return (yp,ys)


from sklearn.metrics import roc_auc_score
from sklearn.metrics import r2_score
import numpy as np


def rmse(labels, predictions):
    return np.sqrt(np.power(labels - predictions, 2).mean())


def mae(labels, predictions):
    return np.abs(labels - predictions).mean()


def r_squared(labels, predictions):
    return r2_score(labels, predictions)


def rocauc(labels, predictions):
    return roc_auc_score(labels, predictions)
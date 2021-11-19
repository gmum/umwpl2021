import numpy as np
import pandas as pd
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem



def get_random_fold_indices(data, n_folds=10):
    data_size = len(data)
    indices = np.random.permutation(data_size)
    fold_size = data_size / n_folds
    return [indices[int(fold_size * i):int(fold_size * (i + 1))] for i in range(n_folds)]


def save_fold_indices(path: str, fold_indices: List[np.ndarray]) -> None:
    np.savez(path, **{f'fold{i}': indices for i, indices in enumerate(fold_indices)})


def load_fold_indices(path: str) -> List[np.ndarray]:
    fold_indices = np.load(path)
    n_folds = len(fold_indices)
    return [fold_indices[f'fold{i}'] for i in range(n_folds)]


def split_data(data: pd.DataFrame, fold_indices: List[np.ndarray], test_fold_idx: int = 0, valid_fold_idx: int = 1) -> Tuple[pd.DataFrame, ...]:
    n_folds = len(fold_indices)
    train_folds = np.array([i for i in range(n_folds) if i not in (test_fold_idx, valid_fold_idx)])
    train_indices = np.concatenate([fold_indices[fold_idx] for fold_idx in train_folds])
    test_indices = fold_indices[test_fold_idx]
    if valid_fold_idx is not None:
        valid_indices = fold_indices[valid_fold_idx]
        return data.iloc[train_indices], data.iloc[valid_indices], data.iloc[test_indices]
    else:
        return data.iloc[train_indices], data.iloc[test_indices]


def cross_validate(data, fold_indices, preprocessing_fn=None, y_column=None):
    n_folds = len(fold_indices)
    for i in range(n_folds):
        train_data, valid_data, test_data = split_data(
            data, 
            fold_indices,
            test_fold_idx=i,
            valid_fold_idx=(i + 1) % n_folds
        )
        if preprocessing_fn is not None:
            train_data, valid_data, test_data = (
                preprocessing_fn(train_data),
                preprocessing_fn(valid_data),
                preprocessing_fn(test_data)
            )
        yield train_data, valid_data, test_data


def load_esol(split_path='../data/esol/split.npz'):
    df = pd.read_csv('https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv')
    fold_indices = load_fold_indices(split_path)
    return df, fold_indices


class Featurizer:
    def __init__(self, y_column, **kwargs):
        self.y_column = y_column
        self.__dict__.update(kwargs)
    
    def __call__(self, df):
        raise NotImplementedError()


class CheatingFeaturizer(Featurizer):
    def __call__(self, df):
        predictions = []
        labels = []
        predictions_column_name = 'ESOL predicted log solubility in mols per litre'
        for i, row in df.iterrows():
            predictions.append(row[predictions_column_name])
            labels.append(row[self.y_column])
        predictions = np.array(predictions)
        labels = np.array(labels)
        return predictions, labels
        

class ECFPFeaturizer(Featurizer):
    def __init__(self, y_column, radius=2, length=1024, **kwargs):
        self.radius = radius
        self.length = length
        super().__init__(y_column, **kwargs)
    
    def __call__(self, df):
        fingerprints = []
        labels = []
        for i, row in df.iterrows():
            y = row[self.y_column]
            smiles = row.smiles
            mol = Chem.MolFromSmiles(smiles)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=self.length)
            fingerprints.append(fp)
            labels.append(y)
        fingerprints = np.array(fingerprints)
        labels = np.array(labels)
        return fingerprints, labels
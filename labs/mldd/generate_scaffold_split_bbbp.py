from collections import defaultdict
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import Chem
import numpy as np
import pandas as pd


def generate_scaffold(smiles, include_chirality=False):
    """return scaffold string of target molecule"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
    else:
        scaffold = ''
    return scaffold

df = pd.read_csv('https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv')
rng = np.random.RandomState(123)
smiles_list = list(df.smiles)
include_chirality = True
frac_valid = 0.1
frac_test = 0.1

scaffolds = defaultdict(list)
for ind, smiles in enumerate(smiles_list):
    scaffold = generate_scaffold(smiles, include_chirality)
    scaffolds[scaffold].append(ind)

scaffold_sets = rng.permutation(list(scaffolds.values()))

n_total_valid = int(np.floor(frac_valid * len(df)))
n_total_test = int(np.floor(frac_test * len(df)))

train_index = []
valid_index = []
test_index = []

for scaffold_set in scaffold_sets:
    if len(valid_index) + len(scaffold_set) <= n_total_valid:
        valid_index.extend(scaffold_set)
    elif len(test_index) + len(scaffold_set) <= n_total_test:
        test_index.extend(scaffold_set)
    else:
        train_index.extend(scaffold_set)

path = '../data/bbbp/split.npz'
np.savez(path, train=np.array(train_index), valid=np.array(valid_index), test=np.array(test_index))
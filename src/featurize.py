#!/usr/bin/env python
"""
Featurization
"""

import sys
import numpy as np
from deepchem import deepchem
from tdc.benchmark_group import admet_group


# get task
group = admet_group(path='data/')
benchmark = group.get(sys.argv[1])
name = benchmark['name']

 # split data
train_val, test = benchmark['train_val'], benchmark['test']

# featurize
maccskeys = deepchem.feat.MACCSKeysFingerprint()
circular = deepchem.feat.CircularFingerprint()
mol2vec = deepchem.feat.Mol2VecFingerprint()
mordred = deepchem.feat.MordredDescriptors(ignore_3D=True)
rdkit = deepchem.feat.RDKitDescriptors()
pubchem = deepchem.feat.PubChemFingerprint()

maccskeys_train_val = maccskeys.featurize(train_val.iloc[:, 1].to_list())
maccskeys_test = maccskeys.featurize(test.iloc[:, 1].to_list())
circular_train_val = circular.featurize(train_val.iloc[:, 1].to_list())
circular_test = circular.featurize(test.iloc[:, 1].to_list())
mol2vec_train_val = mol2vec.featurize(train_val.iloc[:, 1].to_list())
mol2vec_test = mol2vec.featurize(test.iloc[:, 1].to_list())
mordred_train_val = mordred.featurize(train_val.iloc[:, 1].to_list())
mordred_test = mordred.featurize(test.iloc[:, 1].to_list())
rdkit_train_val = rdkit.featurize(train_val.iloc[:, 1].to_list())
rdkit_test = rdkit.featurize(test.iloc[:, 1].to_list())
pubchem_train_val = pubchem.featurize(train_val.iloc[:, 1].to_list())
pubchem_test = pubchem.featurize(test.iloc[:, 1].to_list())

# process pubchem empty lists
pubchem_train_val_tmp = []
for i in range(len(pubchem_train_val)):
    if len(pubchem_train_val[i]) == 0:
        pubchem_train_val[i] = np.array([0] * 881)
    pubchem_train_val_tmp.append(pubchem_train_val[i].tolist())

pubchem_test_tmp = []
for i in range(len(pubchem_test)):
    if len(pubchem_test[i]) == 0:
        pubchem_test[i] = np.array([0] * 881)
    pubchem_test_tmp.append(pubchem_test[i].tolist())

pubchem_test = np.array(pubchem_test_tmp)

# combine features
fp_train_val = np.concatenate(
    (
        maccskeys_train_val, circular_train_val, mol2vec_train_val,
        rdkit_train_val, mordred_train_val, pubchem_train_val
    ), axis=1
)

fp_test = np.concatenate(
    (
        maccskeys_test, circular_test, mol2vec_test,
        rdkit_test, mordred_test, pubchem_test
    ), axis=1
)

# convert nan to 0
fp_train_val = np.nan_to_num(fp_train_val, nan=0, posinf=0)
fp_test = np.nan_to_num(fp_test, nan=0, posinf=0)

# save to npy
np.save(open("./features/" + name + "_train_val.npy", "wb"), fp_train_val)
np.save(open("./features/" + name + "_test.npy", "wb"), fp_test)

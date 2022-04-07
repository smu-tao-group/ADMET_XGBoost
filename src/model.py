#!/usr/bin/env python
"""
Model training
"""

import sys
import numpy as np
import xgboost
from sklearn.model_selection import RandomizedSearchCV
from tdc.benchmark_group import admet_group
from config import config


# get task
group = admet_group(path='data/')
benchmark = group.get(sys.argv[1])
name = benchmark['name']
 # split data
train_val, test = benchmark['train_val'], benchmark['test']
y_train_val = train_val.iloc[:, 2].tolist()

fp_train_val = np.load(open("./features/" + name + "_train_val.npy", "rb"))
fp_test = np.load(open("./features/" + name + "_test.npy", "rb"))

# finetune with grid search
xgb = xgboost.XGBRegressor(tree_method='gpu_hist')

params = {
    "n_estimators": [50, 100, 200, 500, 1000],
    "max_depth": [3, 4, 5, 6, 7],
    "learning_rate": [0.01, 0.05, 0.1, 0.2, 0.3],
    "subsample": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    "colsample_bytree": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    "reg_alpha": [0, 0.1, 1, 5, 10],
    "reg_lambda": [0, 0.1, 1, 5, 10]
}

clf = RandomizedSearchCV(
    estimator=xgb, param_distributions=params,
    scoring=config[name], n_iter=100
)
clf.fit(fp_train_val, y_train_val)

# print best parameters
best_params = clf.best_params_
print(best_params)

predictions_list = []
feature_imp_list = []

for random_state in range(5):
    xgb = xgboost.XGBRegressor(
        tree_method='gpu_hist',
        **best_params,
        random_state=random_state
    )
    xgb.fit(fp_train_val, y_train_val)
    pred_xgb = xgb.predict(fp_test)
    # add to predicitons dict
    predictions = {}
    predictions[name] = pred_xgb
    predictions_list.append(predictions)
    # get feature importance
    feature_imp_list.append(xgb.feature_importances_)
    del xgb

print(group.evaluate_many(predictions_list))

# feature importance by fingerprints and descriptors
# ending 0 is added as placeholder
feature_imp = []
feature_size = [167, 2048, 300, 1613, 208, 0]

for feature in feature_imp_list:
    feature_imp_cur = []
    running_size = 0
    for size in feature_size:
        feature_imp_cur.append(
            np.sum(feature[running_size: running_size + size])
        )
        running_size += size
    feature_imp.append(feature_imp_cur)

feature_imp_mean = np.mean(feature_imp, axis=0)

print(
    f"maccskeys: {feature_imp_mean[0]*100:.2f}% ",
    f"circular: {feature_imp_mean[1]*100:.2f}% ",
    f"mol2vec: {feature_imp_mean[2]*100:.2f}% ",
    f"mordred: {feature_imp_mean[3]*100:.2f}% ",
    f"rdkit: {feature_imp_mean[4]*100:.2f}%"
)

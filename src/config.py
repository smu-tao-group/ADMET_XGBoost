#!/usr/bin/env python
"""
Metrics for ADMET tasks
"""

from scipy.stats import spearmanr
from sklearn.metrics import make_scorer, average_precision_score


def spearman_loss_func(y_true, y_pred):
    """spearman metric
    """
    return spearmanr(y_true, y_pred)[0]


spearman = make_scorer(spearman_loss_func)
auprc = make_scorer(average_precision_score)

config = {
    'caco2_wang': 'neg_mean_absolute_error',
    'bioavailability_ma': 'roc_auc',
    'lipophilicity_astrazeneca': 'neg_mean_absolute_error',
    'solubility_aqsoldb': 'neg_mean_absolute_error',
    'hia_hou': 'roc_auc',
    'pgp_broccatelli': 'roc_auc',
    'bbb_martins': 'roc_auc',
    'ppbr_az': 'neg_mean_absolute_error',
    'vdss_lombardo': spearman,
    'cyp2d6_veith': auprc,
    'cyp3a4_veith': auprc,
    'cyp2c9_veith': auprc,
    'cyp2c9_substrate_carbonmangels': auprc,
    'cyp2d6_substrate_carbonmangels': auprc,
    'cyp3a4_substrate_carbonmangels': 'roc_auc',
    'half_life_obach': spearman,
    'clearance_microsome_az': spearman,
    'clearance_hepatocyte_az': spearman,
    'ld50_zhu': 'neg_mean_absolute_error',
    'herg': 'roc_auc',
    'ames': 'roc_auc',
    'dili': 'roc_auc'
}

best_params = {
    'caco2_wang': {
        'subsample': 0.5, 'reg_lambda': 0.1, 'reg_alpha': 0.1,
        'n_estimators': 1000, 'min_child_weight': 1, 'max_depth': 7,
        'learning_rate': 0.05, 'colsample_bytree': 0.9
    },
    'hia_hou': {
        'subsample': 0.7, 'reg_lambda': 5, 'reg_alpha': 0,
        'n_estimators': 200, 'min_child_weight': 1, 'max_depth': 5,
        'learning_rate': 0.01, 'colsample_bytree': 0.9
    },
    'pgp_broccatelli': {
        'subsample': 0.5, 'reg_lambda': 5, 'reg_alpha': 0.1, 
        'n_estimators': 1000, 'min_child_weight': 3, 'max_depth': 4, 
        'learning_rate': 0.01, 'colsample_bytree': 0.7
    },
    'bioavailability_ma': {
        'subsample': 0.7, 'reg_lambda': 0, 'reg_alpha': 1,
        'n_estimators': 500, 'min_child_weight': 3, 'max_depth': 5,
        'learning_rate': 0.01, 'colsample_bytree': 0.6
    },
    'lipophilicity_astrazeneca': {
        'subsample': 0.8, 'reg_lambda': 5, 'reg_alpha': 0.1,
        'n_estimators': 500, 'min_child_weight': 1, 'max_depth': 4,
        'learning_rate': 0.05, 'colsample_bytree': 0.5
    },
    'solubility_aqsoldb': {
        'subsample': 0.6, 'reg_lambda': 1, 'reg_alpha': 10,
        'n_estimators': 500, 'min_child_weight': 1, 'max_depth': 5,
        'learning_rate': 0.05, 'colsample_bytree': 1.0
    },
    'bbb_martins': {
        'subsample': 0.9, 'reg_lambda': 1, 'reg_alpha': 0.1,
        'n_estimators': 1000, 'min_child_weight': 3, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.9
    },
    'ppbr_az': {
        'subsample': 0.5, 'reg_lambda': 5, 'reg_alpha': 0.1,
        'n_estimators': 1000, 'min_child_weight': 5, 'max_depth': 3,
        'learning_rate': 0.05, 'colsample_bytree': 1.0
    },
    'vdss_lombardo': {
        'subsample': 0.5, 'reg_lambda': 10, 'reg_alpha': 10,
        'n_estimators': 200, 'min_child_weight': 1, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.5
    },
    'cyp2c9_veith': {
        'subsample': 1.0, 'reg_lambda': 5, 'reg_alpha': 0,
        'n_estimators': 500, 'min_child_weight': 5, 'max_depth': 6,
        'learning_rate': 0.05, 'colsample_bytree': 0.8
    },
    'cyp2d6_veith': {
        'subsample': 0.9, 'reg_lambda': 0.1, 'reg_alpha': 0.1,
        'n_estimators': 1000, 'min_child_weight': 3, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.8
    },
    'cyp3a4_veith': {
        'subsample': 0.6, 'reg_lambda': 0, 'reg_alpha': 1,
        'n_estimators': 1000, 'min_child_weight': 1, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.5
    },
    'cyp2c9_substrate_carbonmangels': {
        'subsample': 0.6, 'reg_lambda': 5, 'reg_alpha': 0.1,
        'n_estimators': 200, 'min_child_weight': 1, 'max_depth': 3,
        'learning_rate': 0.01, 'colsample_bytree': 1.0
    },
    'cyp2d6_substrate_carbonmangels': {
        'subsample': 0.5, 'reg_lambda': 5, 'reg_alpha': 0.1,
        'n_estimators': 500, 'min_child_weight': 5, 'max_depth': 4,
        'learning_rate': 0.05, 'colsample_bytree': 0.7
    },
    'cyp3a4_substrate_carbonmangels': {
        'subsample': 0.9, 'reg_lambda': 0, 'reg_alpha': 0,
        'n_estimators': 500, 'min_child_weight': 1, 'max_depth': 5,
        'learning_rate': 0.01, 'colsample_bytree': 0.6
    },
    'half_life_obach': {
        'subsample': 0.7, 'reg_lambda': 10, 'reg_alpha': 10,
        'n_estimators': 100, 'min_child_weight': 1, 'max_depth': 7,
        'learning_rate': 0.05, 'colsample_bytree': 0.6
    },
    'clearance_hepatocyte_az': {
        'subsample': 1.0, 'reg_lambda': 10, 'reg_alpha': 1,
        'n_estimators': 1000, 'min_child_weight': 1, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.7
    },
    'clearance_microsome_az': {
        'subsample': 0.9, 'reg_lambda': 10, 'reg_alpha': 1,
        'n_estimators': 500, 'min_child_weight': 1, 'max_depth': 4,
        'learning_rate': 0.01, 'colsample_bytree': 0.6
    },
    'ld50_zhu': {
        'subsample': 0.5, 'reg_lambda': 0.1, 'reg_alpha': 0.1,
        'n_estimators': 1000, 'min_child_weight': 1, 'max_depth': 7,
        'learning_rate': 0.05, 'colsample_bytree': 0.9
    },
    'herg': {
        'subsample': 0.5, 'reg_lambda': 0.1, 'reg_alpha': 5,
        'n_estimators': 500, 'min_child_weight': 3, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.9
    },
    'ames': {
        'subsample': 0.5, 'reg_lambda': 10, 'reg_alpha': 0.1,
        'n_estimators': 500, 'min_child_weight': 1, 'max_depth': 6,
        'learning_rate': 0.01, 'colsample_bytree': 0.9
    },
    'dili': {
        'subsample': 0.5, 'reg_lambda': 10, 'reg_alpha': 0,
        'n_estimators': 1000, 'min_child_weight': 1, 'max_depth': 6,
        'learning_rate': 0.1, 'colsample_bytree': 0.9
    }
}

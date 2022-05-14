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

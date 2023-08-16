import numpy as np
import pandas as pd
import sklearn

from . import types

def standardized_mutual_information(y_true: types.Array, p_label: types.Array) -> float:
    """MI(ground truth, predicted binary labels) / MI(ground truth, ground truth)."""
    return sklearn.metrics.mutual_info_score(y_true, p_label)/sklearn.metrics.mutual_info_score(y_true, y_true)

def evaluate_classification(y_true: types.Array, p_hat: types.Array,
                            p_label: types.Array) -> types.NumberDict:
    """Binary classification metrics."""

    # subset on cells that are not nan (nan are cells that were not assigned to any nhood in Milo)
    milo_not_nan = ~np.isnan(p_hat)
    y_true = y_true[milo_not_nan]
    p_hat = p_hat[milo_not_nan]
    p_label = p_label[milo_not_nan]
    
    return {
        'AUROC': sklearn.metrics.roc_auc_score(y_true, p_hat),
        'AP': sklearn.metrics.average_precision_score(y_true, p_hat),
        'Recall': sklearn.metrics.recall_score(y_true, p_label),
        'Precision': sklearn.metrics.precision_score(y_true, p_label),
        'F1': sklearn.metrics.f1_score(y_true, p_label),
        'MCC' : sklearn.metrics.matthews_corrcoef(y_true, p_label),
        'Standardized MI': standardized_mutual_information(y_true, p_label)
    }

def _helper_metrics(t_p: int, f_p: int, f_n: int) -> types.NumberDict:
    """A bunch of classification metrics."""
    if (t_p + f_p) == 0:
        precision = np.nan
        fdr = np.nan
    else:
        precision = t_p / (t_p + f_p)
        fdr = f_p / (f_p + t_p)

    recall = t_p / (t_p + f_n)
    if (precision + recall) == 0 or (precision + recall) == np.nan:
        f1_score = np.nan
    else:
        f1_score = 2 * precision * recall / (precision + recall)

        
    return {'Precision': precision, 'Recall': recall, 'FDR': fdr, 'F1': f1_score}


def evaluate_de_genes(genes: types.StrArrayDict,
                      true_a: str, true_b: str,
                      pred_a: str, pred_b: str) -> types.NumberDict:
    """"Classification metrics for differential expression genes."""
    all_true = set(genes[true_a]).union(genes[true_b])
    all_pred = set(genes[pred_a]).union(genes[pred_b])
    tp = len(all_true.intersection(all_pred))
    fn = len(all_true.difference(all_pred))
    fp = len(all_pred.difference(all_true))
    info = _helper_metrics(tp, fp, fn)
    info['n found (0)']= len(genes[pred_a])
    info['n found (1)']= len(genes[pred_b])
    return info

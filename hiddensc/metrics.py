import numpy as np
import sklearn

from . import types


def evaluate_classification(y_true: types.Array, p_hat: types.Array,
                            p_label: types.Array) -> types.NumberDict:
    """Binary classification metrics."""
    return {
        'AUCROC': sklearn.metrics.roc_auc_score(y_true, p_hat),
        'AP': sklearn.metrics.average_precision_score(y_true, p_hat),
        'Recall': sklearn.metrics.recall_score(y_true, p_label),
        'Precision': sklearn.metrics.precision_score(y_true, p_label),
        'F1': sklearn.metrics.f1_score(y_true, p_label)
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
    return _helper_metrics(tp, fp, fn)

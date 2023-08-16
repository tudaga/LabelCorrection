import collections
import gc
from typing import Tuple, List

import numpy as np
import scanpy as sc
import scipy.stats
import sklearn.linear_model
import sklearn.svm
from sklearn.cluster import KMeans
from tqdm.auto import tqdm

from . import types


def get_pca(adata: types.AnnData, n_comps: int = 50) -> np.ndarray:
    """Basic pca, assumes X is count data."""
    x_scaled = sc.pp.scale(adata.X, max_value=10)
    return sc.tl.pca(x_scaled, svd_solver='arpack', n_comps=n_comps)

def map_to_binary(values, case_control_labels, case_cond=1):
    
    # simple thresholding works, since the values are either very close to -1 or very close to 1
    # return np.interp(values, [np.min(values),np.max(values)], [0,1])>=0.5
    
    # but for consistency,
    # using kmeans with n_clusters=2, same as with HiDDEN p_hat binarization of the cells in the case_cond
    
    kmeans_values_res = KMeans(n_clusters=2, random_state=0).fit(pd.DataFrame(values))
    mean_values_res_kmeans_label0 = np.mean(values[kmeans_values_res.labels_==0]) 
    mean_values_res_kmeans_label1 = np.mean(values[kmeans_values_res.labels_==1])
    zero_lab_has_lower_mean = mean_values_res_kmeans_label0 < mean_values_res_kmeans_label1

    df_values_clust = pd.DataFrame(values)
    df_values_clust['kmeans'] = 0
    df_values_clust['kmeans'][(case_control_labels==case_cond).values] = [1 if x==int(zero_lab_has_lower_mean) else 0 for x in kmeans_values_res.labels_]
    
    return df_values_clust['kmeans'].values

def kmeans_correction(y_true: np.ndarray, y_pred: np.ndarray, case_cond: int, rand_state: int) -> np.ndarray:
    """Correct predicted labels using a k-means strategy, helps account for potential bi-modality."""
    is_case_cond = y_true == case_cond
    assert np.sum(is_case_cond), f'Found 0 examples of case condition {case_cond}!'
    kmeans_case = sklearn.cluster.KMeans(n_clusters=2, n_init='auto', random_state=rand_state)
    kmeans_case.fit(y_pred[is_case_cond].reshape(-1, 1))
    # We need to correct kmeans since the initialization can flip the labels.
    k_labels = kmeans_case.labels_
    mean_p_hat_kmeans_label0 = np.mean(y_pred[is_case_cond][k_labels == 0])
    mean_p_hat_kmeans_label1 = np.mean(y_pred[is_case_cond][k_labels == 1])
    zero_lab_has_lower_mean = mean_p_hat_kmeans_label0 < mean_p_hat_kmeans_label1
    p_hat = np.zeros_like(y_pred)
    p_hat[is_case_cond] = np.array([1 if x == int(zero_lab_has_lower_mean) else 0 for x in k_labels])
    return p_hat


def logistic_regression(x: np.ndarray, y: np.ndarray, rand_state: int) -> np.ndarray:
    """Basic logistic regression."""
    model = sklearn.linear_model.LogisticRegression(random_state=rand_state, penalty=None)
    model.fit(x, y)
    predicted_prob = model.predict_proba(x)
    return predicted_prob[:, 1]


def logistic_predictions(x: np.ndarray, y: np.ndarray, case_cond: int, rand_state: int) -> Tuple[
    np.ndarray, np.ndarray]:
    """Linear strategy for getting probabilities and predictions."""
    y_prob = logistic_regression(x, y, rand_state)
    y_labels = kmeans_correction(y, y_prob, case_cond, rand_state)
    return y_prob, y_labels


def svm_predictions(x: np.ndarray, y: np.ndarray, case_cond: int, rand_state: int) -> Tuple[np.ndarray, np.ndarray]:
    """An SVM strategy for getting probabilities and predictions."""
    model = sklearn.svm.SVC(kernel='linear', probability=True)
    model.fit(x, y)
    y_prob = model.predict_proba(x)[:, 1]
    y_labels = kmeans_correction(y, y_prob, case_cond, rand_state)
    return y_prob, y_labels


def determine_pcs_heuristic_ks(adata: types.AnnData, orig_label: str="batch",
                               min_pcs: int = 2, max_pcs: int = 60,
                               rand_state: int = 0) -> Tuple[List[int], List[float], List[float]]:
    """Heuristic decision criterion for number of PCs.
    
    Choose the NUM_PCS that maximizes the Kolmogorov-Smirnov test statistic
    comparing the sample distribution of p_hat for HiDDEN label 0 and HiDDEN label 1
    within the case condition.
    """
    x_pca = get_pca(adata, max_pcs)
    results = collections.defaultdict(list)
    for num_pcs in tqdm(np.arange(min_pcs, max_pcs + 1, 1)):
        adata_ks = adata.copy()
        y = adata_ks.obs[orig_label].astype('int').values
        adata_ks.obs[orig_label] = y
        x = x_pca[:, :num_pcs]
        p_hat = logistic_regression(x, y, rand_state)
        adata_ks.obs['p_hat'] = p_hat
        new_labels = kmeans_correction(y, p_hat, 1, rand_state)
        adata_ks.obs['new_labels'] = new_labels

        conditions = [(adata_ks.obs[orig_label] == 0),
                      (adata_ks.obs[orig_label] == 1) & (adata_ks.obs['new_labels'] == 0),
                      (adata_ks.obs[orig_label] == 1) & (adata_ks.obs['new_labels'] == 1)]
        values = ['control', 'case_new0', 'case_new1']
        adata_ks.obs['three_labels_batch_newlabels'] = np.select(conditions, values)

        ks_stat, ks_pval = scipy.stats.ks_2samp(
            adata_ks.obs['p_hat'][adata_ks.obs['three_labels_batch_newlabels'].isin(['control'])].values,
            adata_ks.obs['p_hat'][adata_ks.obs['three_labels_batch_newlabels'].isin(['case_new1', 'case_new0'])].values,
            alternative='greater')

        results['num_pcs'].append(num_pcs)
        results['ks'].append(ks_stat)
        results['ks_pval'].append(ks_pval)

        # Adata can be huge, so we make sure to clean it.
        del adata_ks
        gc.collect()

    return results['num_pcs'], results['ks'], results['ks_pval']

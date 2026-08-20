import numpy as np
import scipy.spatial
import math
import pandas as pd
from itertools import combinations
from persim import wasserstein
from persim.landscapes import (
    PersLandscapeApprox,
    average_approx,
    snap_pl,
    plot_landscape,
    plot_landscape_simple
)
import random
from ripser import ripser
from itertools import product

"""
================================================================================
Aggregation test for intensity functions over a collection of bandwidths

This code is adapted from the code of Schrab et al., 2023.
================================================================================
"""

def Aggtest(
    X,
    Y,
    alpha=0.05,
    kernel="gaussian", # Only the Gaussian kernel is defined in this code. If you want a different kernel, you must input it yourself.
    number_bandwidths=10,
    optimal_bandwidths=True,
    weight_function=None,
    Rff_approx=False,
    B=1000,
    merging_function=np.min,
    seed=42,
    return_dictionary=False,
):
    """
    Two sample Aggtest using persistence diagrams.
    
    Given data from one distribution and data from another distribution,
    return 0 if the test fails to reject the null 
    (i.e. data comes from the same distribution), 
    or return 1 if the test rejects the null 
    (i.e. their intensity functions are different).
    
    Parameters
    ----------
    X: list of persistence diagrams
    Y: list of persistence diagrams
    alpha: scalar
        The value of alpha must be between 0 and 1.
    kernel: str
        The value of kernel must be "gaussian".
        Only the Gaussian kernel is defined in this code. If you want a different kernel, you must input it yourself.
    number_bandwidths: int
        The number of bandwidths per kernel to include in the collection.
    weight_function: function
        function to compute the weight of each point in persistence diagram.
    B: int
        The number of Monte Carlo samples of permutations to generate.
    seed: int 
        Random seed used for the randomness of the permutations.
    return_dictionary: bool
        If true, a dictionary is returned containing for step-by-step results.
    Returns
    -------
    output : int
        0 if the Aggtest fails to reject the null
            (i.e. data comes from the same distribution)
        1 if the Aggtest rejects the null
            (i.e. intensity functions are different)
    dictionary: dict
        Returned only if return_dictionary is True.
        Dictionary containing the overall output of the Aggtest.

    """    
    # Assertions
    n = len(X)
    m = len(Y)
    assert 0 < alpha  and alpha < 1
    assert kernel in (
        "gaussian",
    )
    assert number_bandwidths > 1 and type(number_bandwidths) == int
    assert B>0 and type(B)==int

    assert weight_function is not None

    if not optimal_bandwidths:
        # Coordinate-wise median-based bandwidth collections
        def compute_bandwidths(distances, number_bandwidths):
            if np.min(distances) < 10 ** (-1):
                d = np.sort(distances)
                lambda_min = np.maximum(
                    d[int(np.floor(len(d) * 0.05))],
                    10 ** (-1)
                )
            else:
                lambda_min = np.min(distances)
            lambda_min = lambda_min / 2
            lambda_max = np.maximum(
                np.max(distances),
                3 * 10 ** (-1)
            )
            lambda_max = lambda_max * 2
            power = (lambda_max / lambda_min) ** (
                1 / (number_bandwidths - 1)
            )
            bandwidths = np.array([
                power ** i * lambda_min
                for i in range(number_bandwidths)
            ])
            return bandwidths
        max_samples = 100
        all_dists_1 = []
        all_dists_2 = []
        for diagramX in X[:max_samples]:
            for diagramY in Y[:max_samples]:
                # Skip empty diagrams
                if len(diagramX) == 0 or len(diagramY) == 0:
                    continue
                # Coordinate 1 pairwise distances
                d1 = np.abs(
                    diagramX[:, 0][:, None]
                    - diagramY[:, 0][None, :]
                ).reshape(-1)
                # Coordinate 2 pairwise distances
                d2 = np.abs(
                    diagramX[:, 1][:, None]
                    - diagramY[:, 1][None, :]
                ).reshape(-1)

                all_dists_1.append(d1)
                all_dists_2.append(d2)
        distances_1 = np.concatenate(all_dists_1)
        distances_2 = np.concatenate(all_dists_2)
        bandwidths_1 = compute_bandwidths(
            distances_1,
            number_bandwidths
        )
        bandwidths_2 = compute_bandwidths(
            distances_2,
            number_bandwidths
        )
        # Cartesian product:
        # (lambda_1, lambda_2)
        bandwidths = np.array(
            list(product(bandwidths_1, bandwidths_2)),
            dtype=float
        )
        number_bandwidths = len(bandwidths)


    else: #optimal bandwidth collection which ensures the minimax optimality up to an iterated logarithm factor.
        T = math.ceil(math.log2((n + m) / math.log(math.log(n + m))))

        bandwidth_1d = np.array([2**(-i) for i in range(1, T + 1)])

        bandwidths = np.array(list(product(bandwidth_1d, repeat=2)))

        number_bandwidths = len(bandwidths)

    # Setup permutations   #Efficient permutation method (Schrab et al., 2023, p. 51)
    rs = np.random.RandomState(seed)
    idx = rs.rand(B, n+m).argsort(axis=1)  # (B, n+m): rows of permuted indices
    #11
    v11 = np.concatenate((np.ones(n), -np.ones(m)))  # (n+m, )
    V11i = np.tile(v11, (B, 1))  # (B, n+m)
    V11 = np.take_along_axis(V11i, idx, axis=1)  # (B, n+m): permute the entries of the rows
    V11 = np.vstack([V11, v11])    # (B+1)th entry is the original test statistic (no permutation)
    V11 = V11.transpose()  # (n+m, B+1)
    #10
    v10 = np.concatenate((np.ones(n), np.zeros(m)))
    V10i = np.tile(v10, (B, 1))
    V10 = np.take_along_axis(V10i, idx, axis=1)
    V10 = np.vstack([V10, v10])
    V10 = V10.transpose()
    #01
    v01 = np.concatenate((np.zeros(n), -np.ones(m)))
    V01i = np.tile(v01, (B, 1))
    V01 = np.take_along_axis(V01i, idx, axis=1)
    V01= np.vstack([V01, v01])
    V01 = V01.transpose()

        
    # Step 1: compute the statistical matrix M
    N = number_bandwidths
    M = np.zeros((N, B + 1))
    list_pd = X + Y # list concatenation
    for i in range(number_bandwidths):
        bandwidth = bandwidths[i]
        K = Mat_gram(list_pd,bandwidth,weight_function,kernel,seed=seed,Rff_approx=Rff_approx,num_rff=10**4)
        # set diagonal elements to zero
        np.fill_diagonal(K, 0)
        # compute T permuted values
        M[i] = (
            np.sum(V10 * (K @ V10), 0) * (m - n + 1) / (m * n * (n - 1))
            + np.sum(V01 * (K @ V01), 0) * (n - m + 1) / (m * n * (m - 1))
            + np.sum(V11 * (K @ V11), 0) / (m * n)
        )
    M=M.transpose() #(B+1,N), the (B+1)th row of M corresponds the original statistic.

    # Step 2: compute P-value Matrix
    def computing_rank(vector):
        vector = np.asarray(vector)
        L = len(vector)
        return np.sum(vector[:, None] <= vector[None, :], axis=1) / L
    Pvar_M=np.column_stack([computing_rank(M[:, j]) for j in range(M.shape[1])]) #P-value matrix

    # Step 3: Aggregated Statistics
    Agg_Stats=np.apply_along_axis(merging_function,axis=1,arr=Pvar_M)
    # Step 4: output test result
    pval= np.sum(Agg_Stats<=Agg_Stats[-1])/len(Agg_Stats) #Agg_Stats[-1] corresponds to the identity permutation.

    
    # create rejection dictionary 
    reject_dictionary = {}

    reject_dictionary["Bandwidth"] = bandwidths
    reject_dictionary["Stat_Matrix"] = M
    reject_dictionary["P-value_Matrix"] = Pvar_M
    reject_dictionary["AggStats"] = Agg_Stats
    reject_dictionary["p-value"] = pval
    # Aggregated test rejects if pval is less than or equal to alpha
    reject_dictionary["Aggtest reject"] = True if pval <= alpha else False

    if return_dictionary:
        return int(reject_dictionary["Aggtest reject"]), reject_dictionary
    else:
        return int(reject_dictionary["Aggtest reject"])



def gaussian_kernel_func(x, y, bandwidth):
    bandwidth = np.asarray(bandwidth)
    return np.exp(-0.5 * np.sum(((x - y) / bandwidth)**2))

def Linear(diagram_1, diagram_2, vec_weight_1, vec_weight_2,kernel,bandwidth):
    """
    diagram_1, _2: Input diagrams, np.arraies of shape (n,2)
    vec_weight_1 , _2: Lists of length n consisting of weight values of points in diagram_1, _2, respectively.
    kernel : A string representing the kernel function
    bandwidth : Flaot, a bandwidth of kernel function
    
    return : Float, a linear kernel value between diagram_1 and diagram_2
    """
    if kernel=='gaussian':
        kernel_func = gaussian_kernel_func
    
    return_value = 0.0
    num_point_1 = diagram_1.shape[0] 
    num_point_2 = diagram_2.shape[0]
    for i in range(num_point_1):
        for j in range(num_point_2):
            return_value += (vec_weight_1[i] * vec_weight_2[j]
                  * kernel_func(diagram_1[i, :], diagram_2[j, :],bandwidth))
    return return_value

def function_weight(name_weight, arc_c=1.0, arc_p=5.0, lin_el=1.0,poly_order=2):
    if name_weight == "arctan":
        def func_weight(bd):
            return np.maximum(np.arctan(math.pow((bd[1] - bd[0]), arc_p) / arc_c), 0.0)
    elif name_weight == "linear":
        def func_weight(bd):
            return (bd[1]-bd[0])
    elif name_weight == "Poly":
        def func_weight(bd):
            return (bd[1]-bd[0])**poly_order
    else:  # unweighted
        def func_weight(bd):
            return 1.0
    return func_weight

def vector_weight(function_weight,diagram):
        num_point = diagram.shape[0]
        vec = np.empty(num_point)
        for k in range(num_point):
            vec[k] = function_weight(diagram[k, :])
        return vec

def list_vector_weight(list_pd,function_weight):
    list_weight = []
    for i in range(len(list_pd)):
        list_weight.append(vector_weight(function_weight=function_weight,diagram=list_pd[i]))
    return list_weight

def Mat_gram(list_pd,bandwidth,weight_function,kernel,seed,Rff_approx=False,num_rff=10**4):
    """
    list_pd: List of length (n+m) which contains all input persistence diagrams
    kernel : String representing the kernel function
    bandwidth : a list representing the bandwidth of kernel function
    weight_function : Function, the weight_function for computing weight of each point in each diagram.
    Rff_approx : Boolean, 1 -> Random Fourier Feature approximation for computing the Linear-kernel value between two diagrams.
    num_rff : The number of samples for the Rff approximation
    
    return : An (n+m) * (n+m) matrix whose element is the Linear kernel value between the corresponding two diagrams.
    """
    num_pd = len(list_pd)
    weight_values = list_vector_weight(list_pd,function_weight= weight_function)
    random.seed(seed)
    np.random.seed(seed)
    if Rff_approx and kernel=='gaussian':
        mat_rff = np.empty((num_pd, num_rff))
        Z= np.random.multivariate_normal(
                [0.0, 0.0], [[bandwidth[0] ** (-2.0), 0.0], [0.0, bandwidth[1] ** (-2.0)]], num_rff)
        b= np.random.uniform(0.0,2.0*math.pi,num_rff)
        for k in range(num_pd):
            mat_rff[k,:]= np.dot(weight_values[k],np.sqrt(2)*np.cos(np.inner(list_pd[k],Z)+b))
        mat =  np.inner(mat_rff, mat_rff) / num_rff
    else:
        mat = np.empty((num_pd, num_pd))
        for i in range(num_pd):
            for j in range(i + 1):
                mat[i, j] = Linear(
                    list_pd[i], list_pd[j],
                    weight_values[i], weight_values[j],kernel,bandwidth)
                mat[j, i] = mat[i, j]
                
    return mat



"""
================================================================================
Persistence Diagram permutation test

This implementation follows Robinson and Turner, 2017.
================================================================================
"""

def compute_distance_matrix(pd_list):
    """
    Computes symmetric matrix of all pairwise Wasserstein distances
    """
    n = len(pd_list)
    dist_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = wasserstein(pd_list[i], pd_list[j])
            dist_mat[i, j] = d
            dist_mat[j, i] = d
    return dist_mat

def avg_pairwise_distance(dist_mat, idxs):
    """
    Computes average pairwise distance within the subset indexed by idxs
    """
    if len(idxs) < 2:
        return 0.0
    submat = dist_mat[np.ix_(idxs, idxs)]
    tril = submat[np.tril_indices(len(idxs), k=-1)]
    return np.mean(tril)

def test_statistic_indices(dist_mat, idx_group1, idx_group2):
    """
    Robinson & Turner test statistic based on indices and precomputed distance matrix
    """
    return avg_pairwise_distance(dist_mat, idx_group1) + avg_pairwise_distance(dist_mat, idx_group2)

def permutation_test(pd_group1, pd_group2, num_permutations=1000, seed=42):
    """
    Fast permutation test with precomputed Wasserstein distances
    """
    random.seed(seed)
    np.random.seed(seed)

    all_diagrams = pd_group1 + pd_group2
    n1 = len(pd_group1)
    n_total = len(all_diagrams)

    dist_mat = compute_distance_matrix(all_diagrams)

    idxs = np.arange(n_total)
    idx_g1 = np.arange(n1)
    idx_g2 = np.arange(n1, n_total)

    T_obs = test_statistic_indices(dist_mat, idx_g1, idx_g2)

    T_perm = []
    for _ in range(num_permutations):
        np.random.shuffle(idxs)
        new_g1 = idxs[:n1]
        new_g2 = idxs[n1:]
        t = test_statistic_indices(dist_mat, new_g1, new_g2)
        T_perm.append(t)

    T_perm = np.array(T_perm)
    p_value = (np.sum(T_perm <= T_obs) + 1) / (num_permutations + 1)
    return T_obs, T_perm, p_value

""""
================================================================================================================================================================
Persistence Landscape test

The base code was referenced from 
https://persim.scikit-tda.org/en/latest/notebooks/Differentiation%20with%20Persistence%20Landscapes.html#Establish-the-baseline
================================================================================================================================================================
"""

def permutation_pl_test(pl_list1,pl_list2,seed=42,num_perms=1000):
    avg_pl_list1 = average_approx(pl_list1)
    avg_pl_list2 = average_approx(pl_list2)
    [avg_pl_list1_snapped, avg_pl_list2_snapped] = snap_pl([avg_pl_list1, avg_pl_list2])
    true_diff_pl = avg_pl_list1_snapped - avg_pl_list2_snapped
    significance = true_diff_pl.p_norm(p=2)

    comb_pl =pl_list1  + pl_list2
    sig_count = 0
    random.seed(seed)
    for shuffle in range(num_perms):
        A_indices = random.sample(range(len(pl_list1)+len(pl_list2)),len(pl_list1))
        B_indices = [_ for _ in range(len(pl_list1)+len(pl_list2)) if _ not in A_indices]

        A_pl = [comb_pl[i] for i in A_indices]
        B_pl = [comb_pl[j] for j in B_indices]

        A_avg = average_approx(A_pl)
        B_avg = average_approx(B_pl)
        [A_avg_sn, B_avg_sn] = snap_pl([A_avg,B_avg])

        shuff_diff = A_avg_sn - B_avg_sn
        if (shuff_diff.p_norm(p=2) >= significance): sig_count += 1

    pval = (sig_count+1)/(num_perms+1)
    return pval
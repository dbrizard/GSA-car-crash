#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:41:13 2024

@author: dbrizard
"""

import numpy as np
import pythia as pt


def sobol_function(x, a=None, **kwargs) -> np.ndarray:
    """Sobol function.

    Parameters
    ----------
    x : np.ndarray
        Input values.
    a : np.ndarray | None, optional
        Coefficients, by default None

    Returns
    -------
    :
        Sobol function values.

    Raises
    ------
    ValueError
        Wrong dimension for input `x`.
    ValueError
        Wrong shape of Coefficients `a`.
    """
    if not 0 < x.ndim < 3:
        raise ValueError("Wrong ndim: {}".format(x.ndim))
    if x.ndim == 1:
        x.shape = 1, -1
    if a is None:
        a = np.zeros(x.shape[1])
    elif not a.shape == (x.shape[1],):
        raise ValueError("Wrong shape: {}".format(a.shape))
    return np.prod((abs(4.0 * x - 2.0) + a) / (1.0 + a), axis=1).reshape(-1, 1)


def sobol_sc(a: np.ndarray, dim: int = 1, **kwargs):
    """Sobol function Sobol indices.

    Parameters
    ----------
    a : np.ndarray
        Coefficients.
    dim : int, optional
        Parameter dimension, by default 1.

    Returns
    -------
    :
        Sobol indices of Sobol function.
    """
    sobol = {}
    beta = (1.0 + a) ** (-2) / 3
    var = np.prod(1.0 + beta) - 1.0
    sobol_tuples = pt.index.IndexSet(pt.index.tensor_set([1, 1, 1])).sobol_tuples
    for sdx in sobol_tuples:
        sobol[sdx] = 1.0 / var
        for k in sdx:
            sobol[sdx] *= beta[k - 1]
    if dim > 1:
        return np.array([sobol for _ in range(dim)])
    else:
        return sobol





def target_function(x: np.ndarray) -> np.ndarray:
    """Target function.

    Parameters
    ----------
    x : np.ndarray
    """
    return sobol_function(x, a=a)


if __name__=="__main__":
    
    #%% Pythia example
    if True:
        print("Tutorial 03 - 2D approximation with PC")

        # target function definition
        a = np.array([1, 2, 3])
        
        # analytical sobol coefficients
        sobol_dict = sobol_sc(a=a, dim=len(a))[0]
        sobol_coefficients = np.array(list(sobol_dict.values())).reshape(-1, 1)
        
        # setup pc surrogate
        params = [
            pt.parameter.Parameter(name=f"x_{j+1}", domain=[0, 1], distribution="uniform")
            for j in range(a.size)
        ]
        
        max_dim = 11
        # limit total polynomial degree of expansion terms to 10
        indices = pt.index.simplex_set(len(params), max_dim - 1)
        index_set = pt.index.IndexSet(indices)
        print("multiindex information:")
        print(f"    number of indices: {index_set.shape[0]}")
        print(f"    dimension: {index_set.shape[1]}")
        print(f"    number of sobol indices: {len(index_set.sobol_tuples)}")
        
        N = 10_000
        print(f"generate training data ({N})")
        s = pt.sampler.WLSTensorSampler(params, [max_dim - 1] * len(params))
        x_train = s.sample(N)
        w_train = s.weight(x_train)
        y_train = target_function(x_train)
        
        print("compute pc expansion")
        surrogate = pt.chaos.PolynomialChaos(params, index_set, x_train, w_train, y_train)
        
        # test PC approximation
        N = 1000
        print(f"generate test data ({N})")
        s_test = pt.sampler.ParameterSampler(params)
        x_test = s_test.sample(N)
        y_test = target_function(x_test)
        y_approx = surrogate.eval(x_test)
        
        error_L2 = np.sqrt(np.sum((y_test - y_approx) ** 2) / N)
        error_L2_rel = error_L2 / np.sqrt(np.sum((y_test) ** 2) / N)
        error_max = np.max(np.abs(y_test - y_approx))
        error_max_rel = np.max(np.abs(y_test - y_approx) / np.abs(y_test))
        
        print(f"test error L2 (abs/rel): {error_L2:4.2e} / {error_L2_rel:4.2e}")
        print(f"test error max (abs/rel): {error_max:4.2e} / {error_max_rel:4.2e}")
        
        # compare Sobol indices
        print("Comparison of Sobol indices")
        print(f" {'sobol_tuple':<12} {'exact':<8}  {'approx':<8}  {'abs error':<9}")
        print("-" * 44)
        for j, sdx in enumerate(sobol_dict.keys()):
            print(
                f" {str(sdx):<11} ",  # Sobol index subscripts
                f"{sobol_coefficients[j, 0]:<4.2e} ",
                f"{surrogate.sobol[j, 0]:<4.2e} ",
                f"{np.abs(sobol_coefficients[j, 0] - surrogate.sobol[j, 0]):<4.2e}",
            )

# coding=utf-8
"""
Usefull definitions for some functions in project
"""

import numpy as np

__author__ = "Arnaud KOPP"
__copyright__ = "© 2014-2015 KOPP Arnaud All Rights Reserved"
__credits__ = ["KOPP Arnaud"]
__license__ = "GNU GPL V2.0"
__version__ = "1.0"
__maintainer__ = "Arnaud KOPP"
__email__ = "kopp.arnaud@gmail.com"
__status__ = "Production"

def adjustpvalues(pvalues, method='fdr', n=None):
    """
    returns an array of adjusted pvalues
    Reimplementation of p.adjust in the R package.
    For more information, see the documentation of the
    p.adjust method in R.

    :param pvalues: numeric vector of p-values (possibly with 'NA's).
    :param method: correction method. Valid values are: : fdr(BH), bonferroni, holm, hochberg, BY
    :param n: number of comparisons, must be at least 'length(p)'; only set
    this (to non-default) when you know what you are doing
    """

    if n is None:
        n = len(pvalues)

    if method == "fdr":
        method = "BH"

    # optional, remove NA values
    p = np.array(pvalues, dtype=np.float)
    lp = len(p)

    assert n <= lp

    if n <= 1:
        return p

    if method == "bonferroni":
        p0 = n * p
    elif method == "holm":
        i = np.arange(lp)
        o = np.argsort(p)
        ro = np.argsort(o)
        m = np.maximum.accumulate((n - i) * p[o])
        p0 = m[ro]
    elif method == "hochberg":
        i = np.arange(0, lp)[::-1]
        o = np.argsort(1 - p)
        ro = np.argsort(o)
        m = np.minimum.accumulate((n - i) * p[o])
        p0 = m[ro]
    elif method == "BH":
        i = np.arange(1, lp + 1)[::-1]
        o = np.argsort(1 - p)
        ro = np.argsort(o)
        m = np.minimum.accumulate(float(n) / i * p[o])
        p0 = m[ro]
    elif method == "BY":
        i = np.arange(1, lp + 1)[::-1]
        o = np.argsort(1 - p)
        ro = np.argsort(o)
        q = np.sum(1.0 / np.arange(1, n + 1))
        m = np.minimum.accumulate(q * float(n) / i * p[o])
        p0 = m[ro]
    elif method == "none":
        p0 = p

    return np.minimum(p0, np.ones(len(p0)))
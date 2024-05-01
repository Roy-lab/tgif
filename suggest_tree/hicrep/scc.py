# 230915 Erika Lee: stripping down to only necessary functions for our use
#
#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Usage: hicrep.py
# Description: Compute HiCRep reproducibility stratum-corrected correlation score (SCCS).
# Reference: Genome Res. 2017 Nov;27(11):1939-1949. doi: 10.1101/gr.220640.117
# The algorithm first normalizes the input contact matrices by the total
# number of contacts and then for each chromosome: 1) mean-filter the input
# matrices with an input window size; 2) exclude common zero entries in
# the input matrices; 3) compute the SCC score. It doesn't have the
# procedure to bootstrap the window-size parameter
#
# Distributed under terms of the GNU General Public License v3.0.
import os
import numpy as np
import scipy.sparse as sp
import math
import sys
import warnings
from hicrep.utils import (
    varVstran,
    upperDiagCsr 
    )

def sccByDiag(m1: sp.coo_matrix, m2: sp.coo_matrix, nDiags: int):
    """Compute diagonal-wise hicrep SCC score for the two input matrices up to
    nDiags diagonals


    Args:
        m1 (sp.coo_matrix): input contact matrix 1
        m2 (sp.coo_matrix): input contact matrix 2
        nDiags (int): compute SCC scores for diagonals whose index is in the
        range of [1, nDiags)
    Returns: `float` hicrep SCC scores
    """
    # convert each diagonal to one row of a csr_matrix in order to compute
    # diagonal-wise correlation between m1 and m2
    m1D = upperDiagCsr(m1, nDiags)
    m2D = upperDiagCsr(m2, nDiags)
    nSamplesD = (m1D + m2D).getnnz(axis=1)
    rowSumM1D = m1D.sum(axis=1).A1
    rowSumM2D = m2D.sum(axis=1).A1
    # ignore zero-division warnings because the corresponding elements in the
    # output don't contribute to the SCC scores
    with np.errstate(divide='ignore', invalid='ignore'):
        cov = m1D.multiply(m2D).sum(axis=1).A1 - rowSumM1D * rowSumM2D / nSamplesD
        rhoD = cov / np.sqrt(
            (m1D.power(2).sum(axis=1).A1 - np.square(rowSumM1D) / nSamplesD ) *
            (m2D.power(2).sum(axis=1).A1 - np.square(rowSumM2D) / nSamplesD ))
        wsD = nSamplesD * varVstran(nSamplesD)
        # Convert NaN and Inf resulting from div by 0 to zeros.
        # posinf and neginf added to fix behavior seen in 4DN datasets
        # 4DNFIOQLTI9G and DNFIH7MQHOR at 5kb where inf would be reported
        # as an SCC score
        wsNan2Zero = np.nan_to_num(wsD, copy=True, posinf=0.0, neginf=0.0)
        rhoNan2Zero = np.nan_to_num(rhoD, copy=True, posinf=0.0, neginf=0.0)

    return rhoNan2Zero @ wsNan2Zero / wsNan2Zero.sum()


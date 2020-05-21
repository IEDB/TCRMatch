"""
Defines methods used to score two TCR CDR3Beta sequences using the algorithm defined in arXiv:1205.6031v2

Author: Austin Crinklaw <acrinklaw@lji.org>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""


import numpy as np
from math import sqrt

cdef float p_beta = 0.11387
cdef int p_kmin = 1
cdef int p_kmax = 30

def fmatrix_k1(double [:, :] blm_qij):
    cdef float b_sum
    cdef Py_ssize_t i, j, k, l
    cdef Py_ssize_t x_max = 20, y_max = 20
    marg = np.zeros(x_max, dtype=float)
    cdef double[:] marg_view = marg
    k1 = np.zeros((x_max,y_max), dtype=float)

    # Normalize matrix by marginal frequencies
    for i in range(x_max):
        b_sum = 0.0
        for j in range(y_max):
            b_sum += blm_qij[i][j];
        marg_view[i] = b_sum

    # Calculate K1
    for k in range(x_max):
        for l in range(y_max):
            k1[k][l] = blm_qij[k][l] / (marg_view[k] * marg_view[l])
            k1[k][l] = pow(k1[k][l], p_beta)

    return k1

cdef k2_prod(int[:] s1, int[:] s2, int start1, int start2, int k, double[:, :] k1):
    cdef float k2 = 1
    cdef int x = 0
    cdef Py_ssize_t i1, i2

    while x < k:
        i1 = s1[x+start1]
        i2 = s2[x+start2]
        k2 *= k1[i1][i2]
        x += 1

    return k2

def k3_sum_mod(int[:] s1, int[:] s2, double [:, :] k1, int l1, int l2):
    cdef float k3 = 0.0
    cdef int start1, start2
    
    k = p_kmin
    while k <= p_kmax:
        if k == 2:
            if k3 < 170:
                return 0
        start1 = 0
        while start1 <= (l1 - k):
            start2 = 0
            while start2 <= (l2 - k):
                prod = k2_prod(s1, s2, start1, start2, k, k1)
                k3 += prod
                start2 += 1
            start1 += 1
        k += 1
  
    return k3

def k3_sum(int[:] s1, int[:] s2, double [:, :] k1, int l1, int l2):
    cdef float k3 = 0.0
    cdef int start1, start2
    
    k = p_kmin
    while k <= p_kmax:
        start1 = 0
        while start1 <= (l1 - k):
            start2 = 0
            while start2 <= (l2 - k):
                prod = k2_prod(s1, s2, start1, start2, k, k1)
                k3 += prod
                start2 += 1
            start1 += 1
        k += 1
  
    return k3

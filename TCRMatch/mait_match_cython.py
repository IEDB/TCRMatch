"""
Scores two TCR CDR3Beta sequences using the algorithm defined in arXiv:1205.6031v2

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
import re
import sys
from math import sqrt
from TCRMatch.mait_match import *
import multiprocessing as mp
import functools
import os

# Global constants required
p_beta = 0.11387
p_kmin = 1
p_kmax = 30
package_directory = os.path.dirname(os.path.abspath(__file__))
blosum_mat = os.path.join(package_directory, 'data/blosum62.qij')
iedb_data = os.path.join(package_directory, 'data/iedb_tcr.tsv')


class Peptide:
    """Simple peptide class containing some things to make life easier

    Attributes:
        sequence: peptide sequence
        i: array of integers that maps amino acid position to alphabet array
        aff: normalization factor (kernel3 with self against self)
    """
    def __init__(self, seq):
        self.seq = seq
        self.i = np.zeros(len(seq), dtype=np.int32)
        self.aff = 0.0


def invalid_seq(seq):
    """Checks validity of amino acid sequence with regex
    
    Args:
        seq: amino acid string
    """
    pattern = re.compile(r'[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y]')
    if len(seq) > 0 and not pattern.findall(seq):
        return False
    else:
        return True


def read_pep_list(file):
    """Reads in and validates amino acid sequences
    Args:
        file: string representing file name and location
    """
    seqlist = []
    with open(file, "r") as inf:
        for line in inf:
            seqlist.append(line.rstrip())

    peplist = []
    for seq in seqlist:
        if invalid_seq(seq):
            raise ValueError(
                'Invalid sequence: {} contains non-amino acid characters'.
                format(seq))
        peplist.append(Peptide(seq))

    return peplist


def read_blossummat_qij():
    """Read in blosum frequency matrix for use in kernel methods
    """
    with open(blosum_mat, "r") as inf:
        linelist = inf.readlines()

    mat = np.zeros((20, 20), dtype=float)
    alphabet = ['' for i in range(20)]

    for line in linelist:
        line = line.rstrip()
        if len(line) == 0: continue
        if line[0] == "#": continue
        if line[:4] == '   A':
            for i in range(20):
                alphabet[i] = line[i * 7 + 3]
            j = 0
        else:
            tvec = line.split(" ")
            for i in range(len(tvec)):
                mat[j][i] = tvec[i]
                mat[i][j] = tvec[i]

            j += 1

    return mat, alphabet


def get_scores(peplist1, peplist2, k1, thresh=0.7):
    """Retrieve kernel 3 scores for two lists of peptides
    
    Args:
        peplist1: list of input peptide sequences
        peplist2: list of IEDB peptide sequences to score against
    """
    scores = []
    for pep1 in peplist1:
        for pep2 in peplist2:
            if pep1.seq == pep2.seq:
                continue
            k3_sco = k3_sum(pep1.i, pep2.i, k1, len(pep1.seq), len(pep2.seq))
            sco = k3_sco / sqrt(pep1.aff * pep2.aff)
            if sco > thresh:
                scores.append((pep1.seq, pep2.seq, sco))
    return scores


def tcrmatch(in_name, out_name, n_threads=1):
    """Primary driver for TCRMatch algorithm

    Args:
        in_name: input filename containing amino acid sequences separated by new lines
        out_name: string where we want to save the file
        n_threads: number of cores on CPU to use
    """
    # Read blosum matrix
    blm_qij, alphabet = read_blossummat_qij()
    k1 = fmatrix_k1(blm_qij)

    # Read list of peptides from files
    peplist1 = read_pep_list(in_name)
    peplist2 = read_pep_list(iedb_data)

    # Encode as blosum indices, calculate normalization factor for k3
    for pep in peplist1:
        for j in range(len(pep.seq)):
            ix = alphabet.index(pep.seq[j])
            pep.i[j] = ix
        pep.aff = k3_sum(pep.i, pep.i, k1, len(pep.seq), len(pep.seq))

    for pep in peplist2:
        for j in range(len(pep.seq)):
            ix = alphabet.index(pep.seq[j])
            pep.i[j] = ix
        pep.aff = k3_sum(pep.i, pep.i, k1, len(pep.seq), len(pep.seq))

    if n_threads == 1:
        res = get_scores(peplist1, peplist2, k1)\

        with open(out_name, "w") as outfile:
            outfile.write("input_sequence\tmatch_sequence\tscore\n")
            for tup in res:
                outfile.write("{}\t{}\t{}\n".format(tup[0], tup[1], tup[2]))
    else:
        # Multiprocessing block - use functools partial to use multiparameters
        n_cpus = n_threads
        pool = mp.Pool(processes=n_cpus)
        res = pool.map(functools.partial(get_scores, peplist2=peplist2, k1=k1),
                       np.array_split(peplist1, n_cpus))
        pool.close()
        pool.join()

        with open(out_name, "w") as outfile:
            outfile.write("input_sequence\tmatch_sequence\tscore\n")
            for chunk in res:
                for tup in chunk:
                    outfile.write("{}\t{}\t{}\n".format(
                        tup[0], tup[1], tup[2]))

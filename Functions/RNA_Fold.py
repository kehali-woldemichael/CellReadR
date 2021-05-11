#!/usr/bin/python

import RNA
import pickle
import math
from Bio.Seq import Seq

def compute_mfeFreq(sequence):
    # create a fold_compound object for the current sequence
    fc = RNA.fold_compound(sequence)

    # compute the MFE and corresponding structure
    (mfe_struct, mfe) = fc.mfe()


    # compute partition function
    (bp_propensity, dG) = fc.pf()

    # compute frequency of MFE structure (the 'hard' way)
    kT = RNA.exp_param().kT / 1000.

    prob_mfe = math.exp((dG - mfe) / kT)

    return prob_mfe, mfe

seq = str(pickle.load(open('seqObject.p', 'rb')))

probMfe, mfe = compute_mfeFreq(seq)


# print(str(seq))
print(100*probMfe)





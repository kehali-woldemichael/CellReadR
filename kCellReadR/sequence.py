# For working with sequence objects
from Bio import SeqIO
import re
import pandas as pd
import numpy as np

# Returns GC content
def metric_gcContent(sequence):
    return round((sequence.count("G") + sequence.count("C"))/(len(sequence)),3)

# Checking for in frame TGG and ATG (both number and indices of occurances)
def return_inFrame(sequence, choice):
    # Definnig stop codons
    stopCodons = ['TAG', 'TAA', 'TGA']

    # Generating list of codons in sequence
    strSeq = str(sequence)
    codons = [strSeq for strSeq in re.split(r'(\w{3})', strSeq) if strSeq]

    # Number of in frame TGG and ATG
    num_inF_TGG = codons.count('TGG')
    num_inF_ATG = codons.count('ATG')
    num_inF_Stop = codons.count(stopCodons[0]) + \
                    codons.count(stopCodons[1]) + \
                    codons.count(stopCodons[2])

    # Indices of TGG, ATG, and defined stop codons
    indicesTGG = \
        np.array([key for key, val in enumerate(codons) if val == 'TGG'])*3
    indicesATG = \
        np.array([key for key, val in enumerate(codons) if val == 'ATG'])*3
    indiciesStop = \
        np.array([key for key, val in enumerate(codons) if val in stopCodons])*3

    if choice == 'all':
        return num_inF_TGG, num_inF_ATG, num_inF_Stop, \
            indicesTGG, indicesATG, indiciesStop
    if choice == 'numTGG': return num_inF_TGG

# Return sesRNAs that are in CDS
def check_inCDS(sequence, searchSequence, isoForm, typeSes):
    if typeSes == 'Reverse':
        return 0 != \
            searchSequence[isoForm].seq.count(sequence.reverse_complement())
    elif typeSes == 'Complement':
        return 0 != searchSequence[isoForm].seq.count(sequence.complement())

# Checks if continuous open reading frame by translating to stop ...
def check_cORF(sequence):
    return len(sequence.translate(to_stop=True)) == len(sequence)/3



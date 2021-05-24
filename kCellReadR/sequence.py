
from kCellReadR.paths import *
# For working with sequence objects
import re
import os 
import pandas as pd
import numpy as np

# For working with sequence objects
from Bio import SeqIO
from Bio.Seq import Seq
# For creating SeqRecord objects

from Bio.SeqRecord import SeqRecord
def metric_gcContent(sequence):
    """Returns GC content"""
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

def check_cORF(sequence):
    """Checks if continuous open reading frame by translating to stop ..."""
    return len(sequence.translate(to_stop=True)) == len(sequence)/3

def return_Complements(seqRecords):
    """Simple function to generate reverse complement and complement 
    Given Seqrecord (also takes Seq objecs as well ... ended up a bit convoluted)"""
    if type(seqRecords) == list:
        reverseComplement_seqRecords = []
        complement_seqRecords = []
        for record in seqRecords:
            reverseComplement_seqRecords.append(record.reverse_complement())
            complement_seqRecords.append(SeqRecord(Seq(str(record.reverse_complement().seq)[::-1]), description = 'Complement'))
        return reverseComplement_seqRecords, complement_seqRecords
    else:
        return seqRecords.reverse_complement(), SeqRecord(Seq(str(seqRecords.reverse_complement().seq)[::-1]), description = 'Complement')


def check_inExonVariants(sesRNA, speciesName, geneName, variantTable, seqDirection):
    """Function for checking if sesRNA how many of exons (checks total and only protein coding)"""
    seqDir_Path = ensembl_BasePath + '/' + speciesName.replace(' ', '_') 
    exonPartialPath = seqDir_Path + '/' + geneName + '_exons_'
    
    exonVariant_Count = 0 
    exonVariant_proteinCoding_Count = 0
    exonVariants = []
    i = 0 
    # Sorting so that goes through exons ... 
    for entry in sorted(os.scandir(seqDir_Path), key=lambda e: e.name):
        i += 1
        if exonPartialPath in entry.path:
            exonVariant_Count += 1 
            temp_seqRecord = list(SeqIO.parse(entry.path, "fasta"))
            if len(temp_seqRecord) != 0:
                exonVariants.append(list(SeqIO.parse(entry.path, "fasta")))
                if variantTable['Type'][exonVariant_Count - 1] == 'protein_coding':
                    exonVariant_proteinCoding_Count += 1
                
    inExon = 0
    inExon_proteinCoding = 0
    for i in range(len(exonVariants)):
        for x in range(len(exonVariants[i])):
            if seqDirection == "Complement" and sesRNA.complement() in exonVariants[i][x].seq:
                inExon += 1
                if variantTable['Type'][i] == 'protein_coding':
                    inExon_proteinCoding += 1 
            if seqDirection == "Reverse" and sesRNA.reverse_complement() in exonVariants[i][x].seq:
                inExon += 1
                if variantTable['Type'][i] == 'protein_coding':
                    inExon_proteinCoding += 1 
                
    return (str(inExon) + "/" + str(exonVariant_Count)), (str(inExon_proteinCoding) + "/" + str(exonVariant_proteinCoding_Count))

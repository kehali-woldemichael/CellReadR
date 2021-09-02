
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

    num_inF_TTGG = 0 
    num_inF_TGGA = 0 
    num_inF_TTGGA = 0 
    for i in indicesTGG: 
        if sequence[i-1:i+3] == "TTGG": num_inF_TTGG += 1 
        if sequence[i:i+4] == "TGGA": num_inF_TGGA += 1 
        if sequence[i-1:i+4] == "TTGGA": num_inF_TTGGA += 1 

    if choice == 'all':
        return num_inF_TGG, num_inF_TTGG, num_inF_TGGA, num_inF_TTGGA, num_inF_ATG, num_inF_Stop, \
            indicesTGG, indicesATG, indiciesStop
    if choice == 'numTGG': return num_inF_TGG
    if choice == 'numStop': return num_inF_Stop, indiciesStop

# Return sesRNAs that are in CDS
def check_inSearchSeq(sequence, searchSequence, typeSes):
    if typeSes == 'Reverse':
        return 0 != \
            searchSequence[0].seq.count(sequence.reverse_complement())
    elif typeSes == 'Complement':
        return 0 != searchSequence[0].seq.count(sequence.complement())

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
    cdsPartialPath = seqDir_Path + '/' + geneName + '_cds_'
    cdnaPartialPath = seqDir_Path + '/' + geneName + '_cdna_'
    
    exonVariant_proteinCoding_Count = len(variantTable[variantTable == 'protein_coding'])

    exonVariants = []
    cdsVariants = []
    cdnaVariants = []

    # Sorting so that goes through exons ... 
    for entry in sorted(os.scandir(seqDir_Path), key=lambda e: e.name):
        if exonPartialPath in entry.path:
            exonVariants.append(list(SeqIO.parse(entry.path, "fasta")))
        if cdsPartialPath in entry.path:
            cdsVariants.append(list(SeqIO.parse(entry.path, "fasta")))
        if cdnaPartialPath in entry.path:
            cdnaVariants.append(list(SeqIO.parse(entry.path, "fasta")))
                
    inExon = 0
    inExon_proteinCoding = 0
    inCDS = 0 
    inCDNA = 0 

    # Loop over total number of variants 
    # Check if sesRNA sequence in exon, protein coding exon, and cDNA 
    for i in range(len(variantTable)):
        rC_sesRNA = sesRNA.reverse_complement()

        # Check if sesRNA in exon variant (with seperate count for protein coding exons) 
        for x in range(len(exonVariants[i])):
            if rC_sesRNA in exonVariants[i][x].seq:
                inExon += 1
                if variantTable[i] == 'protein_coding':
                    inExon_proteinCoding += 1 
        if rC_sesRNA in cdnaVariants[i][0].seq: inCDNA += 1

    # Separte loop to check against CDS 
    for i in range(exonVariant_proteinCoding_Count):
        if rC_sesRNA in cdsVariants[i][0].seq: inCDS += 1
                
    return (str(inExon) + "/" + str(len(exonVariants))),(str(inExon_proteinCoding) + "/" + str(exonVariant_proteinCoding_Count)), (str(inCDS) + "/" + str(len(cdsVariants))), (str(inCDNA) + "/" + str(len(cdnaVariants)))

# Given sequence ... converts to in frame TGGs to TAGs and in frame stops so that first 'T' becomes 'G'
# Had to be careful to only work with in frame codons ... initally had made the mistake to just use string.replace ... this would change out of frame codons as well 
def convert_DNA(sequence, numberConvert):
    # Converting to string object for manipulation 
    strSeq = str(sequence)
    # Generating in frame object variables 
    num_inF_TGG, num_inF_TTGG, num_inF_TGGA, num_inF_TTGGA, num_inF_ATG, num_inF_Stop, indices_inF_TGG, \
        indices_inF_ATG, indices_inF_Stop = \
        return_inFrame(strSeq, 'all')
    print(num_inF_TGG)
    # print(num_inF_Stop)

    # Replacing in frame stop codons in sequence 
    for stop in indices_inF_Stop: 
        stopPairs = [("TAG", "GAG"), ("TAA", "GAA"), ("TGA", "GGA")]
        stopSeq = strSeq[stop:stop+3]
        [stopSeq := stopSeq.replace(a, b) for a, b in stopPairs]
        strSeq = strSeq[:stop] + stopSeq + strSeq[stop+3:]
    
    # Setting number convert to all if 'All' selected as number of TGG to convert 
    if numberConvert == 'All': numberConvert = num_inF_TGG
    # Converting TGG's ... up to number set ... and in order from starting with most central 
    # Sorts indicees by distance from center 
    sorted_indices_centralTGG = np.array(sorted(indices_inF_TGG - (len(strSeq)/2), key = abs)) + (len(strSeq)/2)
    # Converts in frame TGG's ... starting from most central TGG ... up to limit set by numberConvert 
    for i in range(numberConvert):
        currentIndex = int(sorted_indices_centralTGG[i])
        strSeq = strSeq[:currentIndex] + 'TAG' + strSeq[currentIndex+3:]
    # Returns RNA 
    return Seq(strSeq).transcribe()

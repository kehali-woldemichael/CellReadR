# For working with sequence objects
from Bio.Seq import Seq
from Bio import SeqIO
# For calling cmd fuctions
import subprocess
# For checiking if reference sequecnes exist
import os
from kCellReadR.sequence import *


# Importing paths
from kCellReadR.paths import *

# Loads reference sequences from bsubsequenceiomaRT Output folder
def load_referenceSequences(geneName, species):
    martBasePath = martBase + species

    # Generates reference sequence files if necessary
    testPath = martBasePath + '/Reverse_' + geneName + '.fasta'
    if os.path.isfile(testPath) == False:
        # Defining shell command
        command = 'Rscript ' + path_RScript + ' ' + geneName + ' ' + species
        # Running shell command
        refSeq = subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    # Loading sequences for reverse complement gene exons
    rC_fileName = martBasePath + '/Reverse_' + geneName + '.fasta'
    rC_exon_records = list(SeqIO.parse(rC_fileName, "fasta"))

    # Loading sequences for complement gene exons
    C_fileName = martBasePath + '/Complement_' + geneName + '.fasta'
    C_exon_records = list(SeqIO.parse(C_fileName, "fasta"))

    # Loading sequences for gene CDS
    CDS_fileName = martBasePath + '/CDS_' + geneName + '.fasta'
    CDS = list(SeqIO.parse(CDS_fileName, "fasta"))

    # Loading sequences for gene CDS
    cDNA_fileName = martBasePath + '/cDNA_' + geneName + '.fasta'
    cDNA = list(SeqIO.parse(cDNA_fileName, "fasta"))

    return rC_exon_records, C_exon_records, CDS, cDNA

def seq_returnEntrez(sequenceID, retType):
    with Entrez.efetch(
        db="nucleotide", rettype=retType, retmode="text", id=sequenceID
    ) as handle:
        seqRecord = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

    handle = Entrez.efetch(db="nucleotide", id=sequenceID, rettype=retType, retmode="text")

    return seqRecord, handle

# Given sequence ... converts to in frame TGGs to TAGs and in frame stops so that first 'T' becomes 'G'
# Had to be careful to only work with in frame codons ... initally had made the mistake to just use string.replace ... this would change out of frame codons as well
def convert_DNA(sequence, numberConvert):
    # Converting to string object for manipulation
    strSeq = str(sequence)
    # Generating in frame object variables
    num_inF_TGG, num_inF_ATG, num_inF_Stop, indicesTGG, indicesATG, indicesStop = return_inFrame(Seq(strSeq), 'all')
    print(num_inF_TGG)
    # print(num_inF_Stop)

    # Replacing in frame stop codons in sequence
    for stop in indicesStop:
        stopPairs = [("TAG", "GAG"), ("TAA", "GAA"), ("TGA", "GGA")]
        stopSeq = strSeq[stop:stop+3]
        [stopSeq := stopSeq.replace(a, b) for a, b in stopPairs]
        strSeq = strSeq[:stop] + stopSeq + strSeq[stop+3:]

    # Setting number convert to all if 'All' selected as number of TGG to convert
    if numberConvert == 'All': numberConvert = num_inF_TGG
    # Converting TGG's ... up to number set ... and in order from starting with most central
    # Sorts indicees by distance from center
    sorted_indices_centralTGG = np.array(sorted(indicesTGG - (len(strSeq)/2), key = abs)) + (len(strSeq)/2)
    # Converts in frame TGG's ... starting from most central TGG ... up to limit set by numberConvert
    for i in range(numberConvert):
        currentIndex = int(sorted_indices_centralTGG[i])
        strSeq = strSeq[:currentIndex] + 'TAG' + strSeq[currentIndex+3:]
    # Returns RNA
    return Seq(strSeq).transcribe()

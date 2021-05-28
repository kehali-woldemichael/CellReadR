# For working with sequence objects
from Bio import SeqIO
from Bio.Seq import Seq
# For creating SeqRecord objects
from Bio.SeqRecord import SeqRecord

# For calling cmd fuctions
import subprocess
# For checiking if reference sequecnes exist
import os
from kCellReadR.sequence import *


# Importing paths
from kCellReadR.paths import *
from kCellReadR.ensembl import *

# Loads reference sequences from bsubsequenceiomaRT Output folder
def load_referenceSequences(species, geneName, spliceVariant):
    """Generates reference sequence files if necessary"""
    
    spliceVariant = str(spliceVariant)
    save_speciesName = species.replace(" ", "_")
    # Base path for sequences for gene 
    gene_BasePath = ensembl_BasePath + '/' + save_speciesName + '/' + geneName 
    # Path at which to check if sequences downloaded alread 
    test_BasePath = gene_BasePath + '_exons_' + spliceVariant + '_' + save_speciesName + '.fasta'
    # Only run if sequences not already downloaded 
    # print('Downloading' + test_BasePath)
    if os.path.isfile(test_BasePath) == False:
        # Downloading sequences 
        download_ensemblSequences(species, geneName)
        # print('Downloaded')

    # Loading exons file 
    exon_fileName = gene_BasePath + '_exons_' + spliceVariant + '_' + save_speciesName + '.fasta'
    exon_records = list(SeqIO.parse(exon_fileName, "fasta"))
    rC_exon_records, C_exon_records = return_Complements(exon_records)

    # Loading sequences for gene CDS
    CDS_fileName = gene_BasePath + '_cds_' + spliceVariant + '_' + save_speciesName + '.fasta'
    CDS = list(SeqIO.parse(CDS_fileName, "fasta"))

    # Loading sequences for gene CDS
    cDNA_fileName = gene_BasePath + '_cdna_' + spliceVariant + '_' + save_speciesName + '.fasta'
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

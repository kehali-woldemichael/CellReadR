# For working with sequence objects
from Bio import SeqIO
from Bio.Seq import Seq
# For creating SeqRecord objects
from Bio.SeqRecord import SeqRecord

# For calling cli fuctions
import subprocess
# For saving sesRNAs 
from datetime import datetime
# For checiking if reference sequecnes exist
import os

# Importing paths
from kCellReadR.sequence import *
from kCellReadR.paths import *
from kCellReadR.ensembl import *

# Loads reference sequences from bsubsequenceiomaRT Output folder
def load_referenceSequences(species, geneName, spliceVariant):
    """Generates reference sequence files if necessary"""
    
    spliceVariant = str(spliceVariant)
    save_speciesName = species.replace(" ", "_")
    # Base path for sequences for gene 
    gene_BasePath = f"{ensembl_BasePath}/{save_speciesName}/{geneName}"
    # Path at which to check if sequences downloaded alread 
    test_BasePath = f"{gene_BasePath}_exons-{spliceVariant}_{save_speciesName}.fasta"
    # Only run if sequences not already downloaded 
    # print('Downloading' + test_BasePath)
    if os.path.isfile(test_BasePath) == False:
        # Downloading sequences 
        download_ensemblSequences(species, geneName)
        # print('Downloaded')

    # Loading exons file 
    exon_fileName = f"{gene_BasePath}_exons-{spliceVariant}_{save_speciesName}.fasta"
    exon_records = list(SeqIO.parse(exon_fileName, "fasta"))
    rC_exon_records, C_exon_records = return_Complements(exon_records)

    # Loading introns file 
    intron_fileName = f"{gene_BasePath}_introns-{spliceVariant}_{save_speciesName}.fasta"
    intron_records = list(SeqIO.parse(intron_fileName, "fasta"))
    rC_intron_records, C_intron_records = return_Complements(intron_records)

    # Loading sequences for gene CDS
    CDS_fileName = f"{gene_BasePath}_cds-{spliceVariant}_{save_speciesName}.fasta"
    CDS = list(SeqIO.parse(CDS_fileName, "fasta"))

    # Loading sequences for gene CDS
    cDNA_fileName = f"{gene_BasePath}_cdna-{spliceVariant}_{save_speciesName}.fasta"
    cDNA = list(SeqIO.parse(cDNA_fileName, "fasta"))

    # Loading sequences for gene genomic 
    genomic_fileName = f"{gene_BasePath}_genomic-{spliceVariant}_{save_speciesName}.fasta"
    genomic = list(SeqIO.parse(genomic_fileName, "fasta"))

    return rC_exon_records, rC_intron_records, CDS, cDNA, genomic
    

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

def save_all_sesRNAs_DNA(sesRNAs, species, gene):
    savePath = basePath + 'Output/sesRNAs/'
    # Generating BioPython directory if does not exist 
    pathlib.Path(savePath).mkdir(parents=True, exist_ok=True)

    # Generate SeqRecord object for each sequence and append to list 
    outputID = f"{gene}_sesRNA_"
    outputDescription = f"sesRNA for {gene}"

    # Generating sequence record objects (for seperate storage)
    outputSeqMulti = []
    n = 1
    for i in sesRNAs:
        outputSeqMulti.append(SeqRecord(i, id = (outputID + str(n)), description = outputDescription))
        n += 1

    # Write output fasta files 
    saveBase = f"{savePath}{species.replace(' ', '_')}_{gene}_sesRNAs_"
    now = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
    saveName = f"{saveBase}{now}.fasta"
    with open(saveName, "w") as output_handle:
        SeqIO.write(outputSeqMulti, output_handle, "fasta")

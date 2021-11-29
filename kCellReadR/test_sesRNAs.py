import pytest
from kCellReadR.sesRNAs import *
from kCellReadR.sequence import *

# @pytest.mark.one 
def test_initial():
    x = 5
    y = 5 
    assert x == y

def test_generate_sesRNAs():
    # Loading exons file 
    exon_fileName = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/kCellReadR/Arc_exons_1_Human.fasta'
    exon_records = list(SeqIO.parse(exon_fileName, "fasta"))
    rC_exon_records, C_exon_records = return_Complements(exon_records)

    sesRNA_length = 192
    targetChoice = 'exon' # exon, cds, genomic
    parameters = parameters_sesRNA(speciesName, geneName,  spliceVariant, 'Reverse', sesRNA_length, 1, 2, 'None', 40, 70, 20, 20)

    target = rC_exon_records 
    test = CDS

    cds_fileName = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/kCellReadR/Arc_cds_1_Human.fasta'

    exonNumber = 1

    variantTable = pd.read_csv('/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/kCellReadR/Arc_Human_variantTable.csv')

    sesSeq, sequenceMetrics, sesRNA_objs = generate_sesRNA(target, test, parameters, exonNumber, variantTable)
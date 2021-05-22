# for generating any necessary directories
import pathlib
import os
import sys
# For working with sequence objects
from Bio.Seq import Seq
from Bio import SeqIO
# For creating SeqRecord objects
from Bio.SeqRecord import SeqRecord

# Importing paths
from kCellReadR.paths import *

def output_temp_sesRNA(sesRNAs_RNA, geneName, sequenceMetrics):
    # Loading RNAfold as RNA
    sys.path.append(pathToFold)
    import _RNA as RNA

    # Writing sequences as seperate fasta files
    i = 1
    all_mfe = []

    for sesRNA in sesRNAs_RNA:
        (ss, mfe) = RNA.fold(str(sesRNAs_RNA[i-1]))
        all_mfe.append(round(mfe, 3))

        # Making sure that single digit number stast with 0 so that files processed in order
        if i < 10: numSes = '0' + str(i)
        else: numSes = str(i)

#        # Defining output name
        outputName = geneName + '_' + numSes
        outputDescription = "sesRNA #" + numSes

        # Original code block w/o any reversing out output fasta for complement
#        outputRecord = SeqRecord(sesRNA, id = outputName, description = outputDescription)
#        outputFull = basePath + '/Output/BioPython/Temp/' + outputName + '.fasta'

        # Writing the sequence reversed if complement
        # Otherwise messes up RNA-RNA interaction calculations
        # Reversing seems to affect mfe ensemble frequency ...
        if sequenceMetrics['TypeSeq'].iloc[i-1] == 'Reverse':
            outputRecord = SeqRecord(sesRNA, id = outputName, description = outputDescription)
            outputFull = basePath + '/Output/BioPython/Temp/' + outputName + '.fasta'
        else:
            # Reversing sequence
            tempSeq = Seq(str(sesRNA[::-1]))
            outputRecord = SeqRecord(tempSeq, id = outputName, description = outputDescription)
            outputFull = basePath + '/Output/BioPython/Temp/' + outputName + '.fasta'

        with open(outputFull, "w") as output_handle:
            SeqIO.write(outputRecord, output_handle, "fasta")

        # Increasing iterator
        i += 1

    return all_mfe

def generate_RNApred(sesRNAs_DNA, sequenceMetrics, geneName):
    # Generating Temp
    tempPath = basePath + 'Output/BioPython/Temp'
    # Creating directory for temp output if does not exist
    pathlib.Path(tempPath).mkdir(parents=True, exist_ok=True)
    # Generating RNA of sesRNA
    sesRNAs_RNA = return_sesRNA_RNA(sesRNAs_DNA)

    # Just making sure to clear Temp folder before starting
    removeCommand = 'rm -rf ' + tempPath + '/*'
    os.system('rm -rf ')

    # Creating temporary fasta files of sesRNAs (RNA)
    all_mfe = output_temp_sesRNA(sesRNAs_RNA, geneName, sequenceMetrics)

    # Add as column in sequence metrics dataframe
    sequenceMetrics['mfe'] = all_mfe
    sequenceMetrics

    return sequenceMetrics

# Converting to RNA for calculating secondary structure
def return_sesRNA_RNA(sesRNAs_DNA):
    sesRNAs_RNA = []
    for i in range(len(sesRNAs_DNA)):
        sesRNAs_RNA.append(sesRNAs_DNA[i].transcribe())
    return sesRNAs_RNA

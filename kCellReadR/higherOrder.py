# for generating any necessary directories
import pathlib
import os
import sys
# For running bash scripts from inside python ... 
import subprocess
# For working with sequence objects
from Bio import SeqIO
from Bio.Seq import Seq
# For creating SeqRecord objects
from Bio.SeqRecord import SeqRecord

# Importing paths
from kCellReadR.paths import *
# For convert DNA
from kCellReadR.sequence import *


def output_temp_sesRNA(sesRNAs_RNA, geneName, sequenceMetrics):
    # Loading RNAfold as RNA
    sys.path.append(pathToFold)
    import _RNA as RNA

    # Writing sequences as seperate fasta files
    all_mfe = []

    for i, sesRNA in enumerate(sesRNAs_RNA, start = 1):
        (ss, mfe) = RNA.fold(str(sesRNAs_RNA[i-1]))
        all_mfe.append(round(mfe, 3))

        # Making sure that single digit number stast with 0 so that files processed in order
        if i < 10: numSes = '0' + str(i)
        else: numSes = str(i)

        # Defining output name
        outputName = geneName + '_' + numSes
        outputDescription = "sesRNA #" + numSes

        # Writing reverse complement
        outputRecord = SeqRecord(sesRNA, id = outputName, description = outputDescription)
        outputFull = basePath + '/Output/BioPython/Temp/' + outputName + '.fasta'

        with open(outputFull, "w") as output_handle:
            SeqIO.write(outputRecord, output_handle, "fasta")

    return all_mfe

def generate_RNApred(sesRNAs_DNA, sequenceMetrics, geneName, numConvert):
    # Generating Temp
    tempPath = basePath + 'Output/BioPython/Temp'
    # Creating directory for temp output if does not exist
    pathlib.Path(tempPath).mkdir(parents=True, exist_ok=True)

    if numConvert == 0:
        # Generating RNA of sesRNA
        sesRNAs_RNA = return_sesRNA_RNA(sesRNAs_DNA)
    elif numConvert > 0: 
        # Version testing with stops converted and central TGG
        sesRNAs_RNA = return_converted_sesRNA_RNA(sesRNAs_DNA, numConvert)

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

# Converting to RNA for calculating secondary structure
def return_converted_sesRNA_RNA(sesRNAs_DNA, numConvert):
    sesRNAs_RNA = []
    for i in range(len(sesRNAs_DNA)):
        sesRNAs_RNA.append(convert_DNA(sesRNAs_DNA[i], numConvert))
    return sesRNAs_RNA


def generate_mfeProb(sequenceMetrics, geneName, species, spliceVariant): 
    """Tests for energetic favorability of secondary structure with RNAfold 
    Args:
        sequenceMetrics ([type]): [description]
        geneName ([type]): [description]
        species ([type]): [description]
        spliceVariant ([type]): [description]

    Returns:
        sequenceMetrics [pandas.DataFrame]: [updated sequence metrics dataframe with mfe frequencies]
    """

    rnaFold_prob = []

    spliceVariant = str(spliceVariant)
    save_speciesName = species.replace(" ", "_")
    gene_BasePath = f"{ensembl_BasePath}/{save_speciesName}/{geneName}"
    CDS_fileName = f"{gene_BasePath}_cds_{spliceVariant}_{save_speciesName}.fasta"
    
    pathTemp = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/Output/BioPython/Temp'
    pathOutTempFold = f"{pathTemp}/temp.out"

    # sorting files in output of scandir 
    for entry in sorted(os.scandir(pathTemp), key=lambda e: e.name):
        # Defining command for RNAfold 
        commandFold = 'RNAfold -p -d2 --noLP < ' + entry.path + ' > ' + pathOutTempFold    
        # Generating RNAfold predictions 
        generateProb = subprocess.run(commandFold, shell=True, stdout=subprocess.PIPE)
        
        # Moving to Temp directory to work on fasta files 
        currentWD = os.getcwd()
        os.chdir('/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/Output/BioPython/Temp')

        # Running script for getting probabilities from RNAfold output file (added to ArchBin btw)
        readProb = subprocess.Popen("rnaFold_prob.sh", shell=True, stdout=subprocess.PIPE)
        returnedProb = readProb.stdout.read()
        # Waiting for last command to finish before storing value in temp.out file 
        readProb.wait()
        # Append frequences ... convert to percentage 
        rnaFold_prob.append(round(float(returnedProb)*100, 3))
        
        # Removing temp.out after finishing each run 
        os.system('rm -rf temp.out')
        os.system('rm -rf temp.csv')
        # Return to initial working directory 
        os.chdir(currentWD)

    # Removing files generated by RNAfold 
    os.system('rm -rf *ss.ps')
    os.system('rm -rf *dp.ps')
    
    # Adding RNA fold mfe ensemble frequency to sequenceMetrics 
    sequenceMetrics['mfeFreq'] = rnaFold_prob
    
    return sequenceMetrics  


def output_RIblast(sequenceMetrics, geneName, testSpecies, spliceVariant, targetName):
    """Runs tests of RNA-RNA interctions with RIblast 

    Args:
        sequenceMetrics ([pandas.DataFrame]): [table with sesRNA metrics]
        geneName ([str]): [Name of gene]
        testSpecies ([str]): [Name of species against which binding is tested]
        spliceVariant ([int]): [Splice variant being tested]
        targetName ([str]): [Type of target to be tested against ... cDNA, CDS, genomic]

    Returns:
        metricsTable_higherOrder [pandas.DataFrame]: [updated metrics table with RNA-RNA binding]
        useful_RIblast [pandas.DataFrame]: [rest of statistics for top two RNA-RNA interactions]
    """

    spliceVariant = str(spliceVariant)
    save_test_speciesName = testSpecies.replace(" ", "_")
    gene_BasePath = ensembl_BasePath + '/' + save_test_speciesName + '/' + geneName 
    
    # Paths to fasta files for test sequence ... depending on given type
    if targetName == 'CDS':
        target_fileName = f"{gene_BasePath}_cds-{spliceVariant}_{save_test_speciesName}.fasta"
    elif targetName == 'cDNA':
        target_fileName = f"{gene_BasePath}_cdna-{spliceVariant}_{save_test_speciesName}.fasta"
    elif targetName == 'genomic':
        target_fileName = f"{gene_BasePath}_genomic-{spliceVariant}_{save_test_speciesName}.fasta"
    print(target_fileName)
    
    # Path to tempRIblast folder 
    path_tempRIblast = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/Output/RIblast/'
    query_Name = f"{path_tempRIblast}{geneName}_db"
    
    # Generating query database 
    commandQuery = f"RIblast db -i {target_fileName} -o {query_Name}"
    os.system(commandQuery)

    # Path and file name for output CSV 
    outputName = f"{path_tempRIblast}{geneName}.csv"
    # Path to directory sesRNA files 
    path_sesRNAs = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/Output/BioPython/Temp'
    
    # Generating pd.DataFrame for storing calculated values 
    columns_RIblast = [' Accessibility Energy', ' Hybridization Energy', ' Interaction Energy', ' BasePair', 
                       ' Accessibility Energy', ' Hybridization Energy', ' Interaction Energy', ' BasePair']
    useful_RIblast =  pd.DataFrame(columns = columns_RIblast)
    
    # Iteratively generating calculations for sesRNA-target interaction 
    # Made sure to go through sesRNA files in order 
    for entry in sorted(os.scandir(path_sesRNAs), key=lambda e: e.name):
        # Running RIblast calculations 
        os.system(f"RIblast ris -i {entry.path} -o {outputName} -d {query_Name}")
        print(entry.path)
        
        # Remove first two lines from CVS to allow for parsing into pandas.Dataframe 
        os.system(f"sed -i 1,2d {outputName}")

        outputCSV =  pd.read_csv(outputName, skiprows=[1])
        # Sorting by hybirzation energy ... have to have extra white space before column name 
        sorted_outputCSV = outputCSV.sort_values(' Hybridization Energy')

        topHybridizationE = sorted_outputCSV[[' Accessibility Energy', ' Hybridization Energy', ' Interaction Energy', ' BasePair']].iloc[0:1]
        secondHybridizationE = sorted_outputCSV[[' Accessibility Energy', ' Hybridization Energy', ' Interaction Energy', ' BasePair']].iloc[1:2]
        temp_RIblast_ouput = pd.concat([topHybridizationE.reset_index(drop=True), secondHybridizationE.reset_index(drop=True)], axis = 1)

        # Appending calculations for current sesRNA values 
        useful_RIblast = useful_RIblast.append(temp_RIblast_ouput)
        
        # Clearing csv 
        os.system(f"rm -rf {outputName}")
        # Clear BioPython temp fasta file for sesRNA 
        os.system(f"rm -rf {entry.path}")
    
    # Clear RIblast Temp directory 
    os.system(f"rm -rf {path_tempRIblast}*")
    

    # Adding RNA-RNA binding statistics for most favorable interaction to rest of sequence metrics table
    metricsTable_higherOrder = pd.concat([sequenceMetrics.reset_index(drop=True), useful_RIblast.reset_index(drop=True).iloc[:, 0:4]], axis = 1)

    return metricsTable_higherOrder, useful_RIblast
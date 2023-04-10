# Paths to be used various functions 

from os.path import basename, dirname, splitext, split
basePath = split(dirname(__file__))[0]

path_outputEnsembl = basePath + '/Output/EnsemblSeq/'
pathToFold = basePath + '/Packages/ViennaRNA_Python3/usr/lib/python3.9/site-packages/RNA'
path_RScript = basePath + '/kCellReadR/refSeq.R'
martBase = basePath + '/Output/biomaRt/'
path_speciesEnsembl = basePath +  '/Reference/SpeciesEnsembl.csv'
ensembl_BasePath = basePath + '/Output/EnsemblSeq'
pathTemp = f'{basePath}/Output/BioPython/Temp'
#path_tempRIblast = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/Output/RIblast/'
path_tempRIblast = f'{basePath}/Output/RIblast/'
# Path to directory sesRNA files 
#path_sesRNAs = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR/Code/Output/BioPython/Temp'
path_sesRNAs = f'{basePath}/Output/BioPython/Temp'

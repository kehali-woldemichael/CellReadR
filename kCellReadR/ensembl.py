# For fetching Ensembl transcript IDs 
import mygene
import requests, sys
import re
import os
import pathlib 
import pandas as pd
# For timing main call 
import time
import sys
import json 
# For sending gene splice variants 
import numpy as np
import jsons

# For converting strings to Seq --> SeqRecord objects and write as fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from Bio import SeqIO


def json_output_spliceVariants():
    '''Ouputs pandas.DataFrame as json ... for node.js Express server
    Splice variant information (bioType)'''

    speciesName = sys.argv[1]
    geneName = sys.argv[2]

    # Generating transcriptIDs for gene
    ensembl_transcriptIDs = return_ensemblTranscriptIDs(speciesName, geneName)
    # Generating transcript splice variant metrics 
    transcriptMetrics = table_transcriptsInfo(ensembl_transcriptIDs)
    # Generating pd.Dataframe
    df = pd.DataFrame(transcriptMetrics)

    # Converting DataFrame to json and dumping it to std.out
    df_json = df.reset_index().to_json(orient="values")
    parsed = jsons.loads(df_json)
    print(jsons.dumps(parsed, indent=4))


def download_ensemblSequences(speciesName = 'noSpecies', geneName = 'noGene'):
    """For getting exons ... for portal
    Given species, gene name, and splice variant number 
    Writes exons, cds, and cDNA of each splice variant separately as single fasta file"""

    # To allow passing parameters either from cli as sys.argv or function paramaters 
    if speciesName == 'noSpecies' and geneName == 'noGene':
        speciesName = sys.argv[1]
        geneName = sys.argv[2]
        # spliceVariant = int(sys.argv[3])

    
    # Get transcript IDs for given gene symbol 
    ensembl_transcriptIDs = return_ensemblTranscriptIDs(speciesName, geneName)

    # Return error if failed to find transript id 
    if ensembl_transcriptIDs == 'Could not find transript ID':
        print('Could not find transcript ID. Manually upload target sequence and submit error.')
    else:
        # id = ensembl_transcriptIDs[spliceVariant]

        save_speciesName = speciesName.replace(" ", "_")
        # Creating directory for ensembl species specific output if does not exist
        # Replaced spaces with _ for naming directories 
        path_species = path_outputEnsembl + "/" + save_speciesName
        pathlib.Path(path_species).mkdir(parents=True, exist_ok=True)

        i = 1

        # Had to split up ... when only one splice variant returns str ... otherwise returns list
        if type(ensembl_transcriptIDs) != str:
            for id in ensembl_transcriptIDs:
                seqNum = str(i)
                outputName_exons = path_species + "/" + geneName + "_exons_" + seqNum + "_" + save_speciesName + ".fasta"
                outputName_introns = path_species + "/" + geneName + "_introns_" + seqNum + "_" + save_speciesName + ".fasta"
                outputName_cds = path_species + "/" + geneName + "_cds_" + seqNum + "_" + save_speciesName + ".fasta"
                outputName_cdna = path_species + "/" + geneName + "_cdna_" + seqNum + "_" + save_speciesName + ".fasta"
                outputName_genomic = path_species + "/" + geneName + "_genomic_" + seqNum + "_" + save_speciesName + ".fasta"

                # Only proceeds if file does not exist
                if os.path.isfile(outputName_exons) != True:

                    # Downloading exon/intron sequences 
                    totalSeq, exons, introns = return_exons(id)
                    # Generating seq record objects 
                    exonRecords = [SeqRecord(Seq(exon), id = geneName, description = "exon") for exon in exons]
                    intronRecords = [SeqRecord(Seq(intron), id = geneName, description = "intron") for intron in introns]
                    # Writing fasta files 
                    SeqIO.write(exonRecords, outputName_exons, "fasta")
                    SeqIO.write(intronRecords, outputName_introns, "fasta")

                    # Writing cds
                    if metric_transcript(id)[2] == 'protein_coding':
                        cds = return_ensemblCDS(id, speciesName)
                        cdsRecords = SeqRecord(Seq(cds), id = geneName, description = "cds")
                        SeqIO.write(cdsRecords, outputName_cds, "fasta")

                    # Writing cdna
                    cdna = return_ensemblCDNA(id, speciesName)
                    cdnaRecords = SeqRecord(Seq(cdna), id = geneName, description = "cds")
                    SeqIO.write(cdnaRecords, outputName_cdna, "fasta")
                    
                    # Writing cdna
                    genomic = return_ensemblGenomic(id, speciesName)
                    genomicRecords = SeqRecord(Seq(genomic), id = geneName, description = "genomic")
                    SeqIO.write(genomicRecords, outputName_genomic, "fasta")

                else:
                    print("Exists")

                
                i += 1
        else:
            outputName_exons = path_species + "/" + geneName + "_exons_" + str(i) + "_" + save_speciesName + ".fasta"
            outputName_introns = path_species + "/" + geneName + "_introns_" + str(i) + "_" + save_speciesName + ".fasta"
            outputName_cds = path_species + "/" + geneName + "_cds_" + str(i) + "_" + save_speciesName + ".fasta"
            outputName_cdna = path_species + "/" + geneName + "_cdna_" + str(i) + "_" + save_speciesName + ".fasta"
            outputName_genomic = path_species + "/" + geneName + "_genomic_" + str(i) + "_" + save_speciesName + ".fasta"

            # Only proceeds if file does not exist
            if os.path.isfile(outputName_exons) != True:

                # Fetching exon/intron records 
                totalSeq, exons, introns = return_exons(ensembl_transcriptIDs)
                exonRecords = [SeqRecord(Seq(exon), id = geneName, description = "exon") for exon in exons]
                intronRecords = [SeqRecord(Seq(intron), id = geneName, description = "intron") for intron in introns]
                SeqIO.write(exonRecords, outputName_exons, "fasta")
                SeqIO.write(intronRecords, outputName_introns, "fasta")

                # Writing cds
                if metric_transcript(ensembl_transcriptIDs)[2] == 'protein_coding':
                    cds = return_ensemblCDS(ensembl_transcriptIDs, speciesName)
                    cdsRecords = SeqRecord(Seq(cds), id = geneName, description = "cds")
                    SeqIO.write(cdsRecords, outputName_cds, "fasta")

                # Writing cdna
                cdna = return_ensemblCDNA(ensembl_transcriptIDs, speciesName)
                cdnaRecords = SeqRecord(Seq(cdna), id = geneName, description = "cds")
                SeqIO.write(cdnaRecords, outputName_cdna, "fasta")

                # Writing genomic
                genomic = return_ensemblGenomic(ensembl_transcriptIDs, speciesName)
                genomicRecords = SeqRecord(Seq(genomic), id = geneName, description = "genomic")
                SeqIO.write(genomicRecords, outputName_genomic, "fasta")
                
            else:
                print("Exists")

    
def metric_transcript(transcriptID):
    """Returns information on transcript"""
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/" + transcriptID + "?expand=1"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if r.json()['biotype'] != 'protein_coding': 
        length = 'no protein'
    else:
        length = r.json()['Translation']['length'] 
 
    return [r.json()['display_name'], r.json()['assembly_name'], r.json()['biotype'], length, str(bool(r.json()['is_canonical']))]

def table_transcriptsInfo(ensembl_transcriptIDs):

    transcriptInfo = [] 

    i = 1
    if type(ensembl_transcriptIDs) != str:
        for transcript in ensembl_transcriptIDs:
            temp_transcriptInfo = []
            temp_transcriptInfo.append(i)
            temp_transcriptInfo.append(transcript)
            temp_transcriptInfo.extend(metric_transcript(transcript))

            transcriptInfo.append(temp_transcriptInfo)
            i += 1

        transcriptMetrics = \
                pd.DataFrame(transcriptInfo, columns = 
                        ['TranscriptNum', 'TranscriptID', 'TranscriptName', 'Assembly', 'Type', 'AA_Length', 'Is_Canonical'])
    else:
        temp_transcriptInfo = []
        temp_transcriptInfo.append(i)
        temp_transcriptInfo.append(ensembl_transcriptIDs)
        temp_transcriptInfo.extend(metric_transcript(ensembl_transcriptIDs))

        transcriptInfo.append(temp_transcriptInfo)
        transcriptMetrics = \
                pd.DataFrame(transcriptInfo, columns = 
                        ['TranscriptNum', 'TranscriptID', 'TranscriptName', 'Assembly', 'Type', 'AA_Length', 'Is_Canonical'])
    

    return transcriptMetrics


def return_ensemblTranscriptIDs(species, geneSymbol):
    """Returns Ensembl stable transcript IDs given species and gene symbol"""

    # Species that can be directly accessed in mygene
    supportedSpecies = ['human', 'mouse', 'rat', 'fruitfly', 'nematode', 'zebrafish', 'thale-cress', 'frog', 'pig']

    # Makes sure that in supported short species names 
    if species.lower() in supportedSpecies:
        mg = mygene.MyGeneInfo()
        symbolSearch = 'symbol:' + geneSymbol 

        tempID_search = mg.query(geneSymbol, scopes='symbol', fields=['ensembltranscript'], species=species)['hits']

        if len(tempID_search) == 0:
            return 'Could not find transript ID'
        else:
            tempID = tempID_search[0]['_id']
            return mg.getgene(tempID, fields = 'ensembl')['ensembl']['transcript']
    # Otherwise much use taxon id
    else:
        scientificSpeciesName, taxonID = return_scientificName(species) 
            
        mg = mygene.MyGeneInfo()
        symbolSearch = 'symbol:' + geneSymbol 

        tempID_search = mg.query(geneSymbol, scopes='symbol', fields=['ensembltranscript'], species=taxonID)['hits']
        if len(tempID_search) == 0:
            return 'Could not find transript ID'
        else:
            tempID = tempID_search[0]['_id']
            return mg.getgene(tempID, fields = 'ensembl')['ensembl']['transcript']



def return_ensemblCDNA(transcriptID, speciesName):
    """Returns sequence of CDS given transript ID"""
    
    server = "https://rest.ensembl.org"
    extBase = "/sequence/id/"
    scientificSpeciesName, taxonID = return_scientificName(speciesName) 
    ext = extBase + transcriptID + "?" + "type=cdna;" + "species=" + scientificSpeciesName 
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    return r.text

def return_ensemblCDS(transcriptID, speciesName):
    """Returns sequence of CDS given transript ID"""
    
    server = "https://rest.ensembl.org"
    extBase = "/sequence/id/"
    scientificSpeciesName, taxonID = return_scientificName(speciesName) 
    ext = extBase + transcriptID + "?" + "type=cds;" + "species=" + scientificSpeciesName 
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    return r.text

def return_ensemblGenomic(transcriptID, speciesName):
    """Returns sequence of CDS given transript ID"""
    
    server = "https://rest.ensembl.org"
    extBase = "/sequence/id/"
    scientificSpeciesName, taxonID = return_scientificName(speciesName) 
    ext = extBase + transcriptID + "?" + "type=genomic;" + "species=" + scientificSpeciesName 
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    return r.text

def split_ExonsIntrons(seq):
    """Generates chunks of upper and lower case sequence letters
    Due to soft mask_feature ... this effectively seperates out exons and intros"""
    return [chunk for chunk in re.split(r"([A-Z]+)", seq) if chunk]


def return_exons(transcriptID):
    """Returns sequences of soft masked total, exons, and introns
    Given Ensembl stable transcript ID"""
    
    server = "https://rest.ensembl.org"
    extBase = "/sequence/id/"
    # Soft mask with type as genomic ... turn introns to lower case 
    ext = extBase + transcriptID + "?" + "type=genomic;mask_feature=soft"

    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    # Splits contiguous upper case and lower case sequence regions (Exons/Upper vs Introns/Lower)
    splitSeq = split_ExonsIntrons(r.text)

    # Seperates out exons and introns from each other 
    exons = [exon for exon in splitSeq if exon.isupper()]
    introns = [intron for intron in splitSeq if intron.islower()]
    
    return r.text, exons, introns


def return_scientificName(speciesName):
    """Returns scientific name in correct format given species common name
    For accessing ensembl rest API"""
    
    # csv file from ensembl of species comman names and etc ... 
    speciesEnsembl = pd.read_csv(path_speciesEnsembl)
    scientificName = speciesEnsembl[speciesEnsembl['Common name'] == speciesName]['Scientific name']
    taxonID = speciesEnsembl[speciesEnsembl['Common name'] == speciesName]['Taxon ID']

    return scientificName.tolist()[0].lower().replace(" ", "_"), taxonID.tolist()[0] 


if __name__ == "__main__":
    from paths import *
    json_output_spliceVariants()
if __name__ != "__main__":
    sys.path.append('/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR_Temp/Code') 
    from kCellReadR.paths import *


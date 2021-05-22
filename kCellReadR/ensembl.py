# For fetching Ensembl transcript IDs 
import mygene
import requests, sys
import re
import os
import pathlib 
import pandas as pd
# For timing main call 
import time
# For sending gene splice variants 
import numpy as np
import jsons

# For converting strings to Seq --> SeqRecord objects and write as fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from Bio import SeqIO


def main():
    """For getting exons ... for portal
    Given species, gene name, and splice variant number 
    Writes exons, cds, and cDNA of each splice variant separately as single fasta file"""
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

        # Creating directory for ensembl species specific output if does not exist
        # Replaced spaces with _ for naming directories 
        path_species = path_outputEnsembl + "/" + speciesName.replace(" ", "_")
        pathlib.Path(path_species).mkdir(parents=True, exist_ok=True)

        i = 1

        # Had to split up ... when only one splice variant returns str ... otherwise returns list
        if type(ensembl_transcriptIDs) != str:
            for id in ensembl_transcriptIDs:
                outputName_exons = path_species + "/" + geneName + "_exons_" + str(i) + "_" + speciesName + "_" + id + ".fasta"
                outputName_cds = path_species + "/" + geneName + "_cds_" + str(i) + "_" + speciesName + "_" + id + ".fasta"
                outputName_cdna = path_species + "/" + geneName + "_cdna" + str(i) + "_" + speciesName + "_" + id + ".fasta"

                # Only proceeds if file does not exist
                if os.path.isfile(outputName_exons) != True:

                    # Writing exons ... only if length greater than 200 bp 
                    totalSeq, exons, introns = return_exons(id)
                    exonRecords = [SeqRecord(Seq(exon), id = geneName, description = "exon") for exon in exons if len(exon) > 200]
                    SeqIO.write(exonRecords, outputName_exons, "fasta")

                    # Writing cds
                    if metric_transcript_bioType(id) == 'protein_coding':
                        cds = return_ensemblCDS(id, speciesName)
                        cdsRecords = SeqRecord(Seq(cds), id = geneName, description = "cds")
                        SeqIO.write(cdsRecords, outputName_cds, "fasta")

                    # Writing cdna
                    cdna = return_ensemblCDNA(id, speciesName)
                    cdnaRecords = SeqRecord(Seq(cdna), id = geneName, description = "cds")
                    SeqIO.write(cdnaRecords, outputName_cdna, "fasta")
                    

                else:
                    print("Exists")

                i += 1
        else:
            outputName_exons = path_species + "/" + geneName + "_exons_" + str(i) + "_" + speciesName + "_" + ensembl_transcriptIDs + ".fasta"
            outputName_cds = path_species + "/" + geneName + "_cds_" + str(i) + "_" + speciesName + "_" + ensembl_transcriptIDs + ".fasta"
            outputName_cdna = path_species + "/" + geneName + "_cdna" + str(i) + "_" + speciesName + "_" + ensembl_transcriptIDs + ".fasta"

            # Only proceeds if file does not exist
            if os.path.isfile(outputName_exons) != True:
                totalSeq, exons, introns = return_exons(ensembl_transcriptIDs)

                # Writing exons ... only if length greater than 200 bp 
                exonRecords = [SeqRecord(Seq(exon), id = geneName, description = "exon") for exon in exons if len(exon) > 200]
                SeqIO.write(exonRecords, outputName_exons, "fasta")

                # Writing cds
                if metric_transcript_bioType(ensembl_transcriptIDs) == 'protein_coding':
                    cds = return_ensemblCDS(ensembl_transcriptIDs, speciesName)
                    cdsRecords = SeqRecord(Seq(cds), id = geneName, description = "cds")
                    SeqIO.write(cdsRecords, outputName_cds, "fasta")

                # Writing cdna
                cdna = return_ensemblCDNA(ensembl_transcriptIDs, speciesName)
                cdnaRecords = SeqRecord(Seq(cdna), id = geneName, description = "cds")
                SeqIO.write(cdnaRecords, outputName_cdna, "fasta")

                
            else:
                print("Exists")

    
def metric_transcript_bioType(transcriptID):
    """Return biotype of specific splice variant"""
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/" + transcriptID + "?expand=1"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    return r.json()['biotype']

def table_transcripsInfo(ensembl_transcriptIDs):
    transcriptNum = []
    transcriptIDs = []
    transcriptType = []
    
    i = 1
    if type(ensembl_transcriptIDs) != str:
        for transcript in ensembl_transcriptIDs:
            transcriptNum.append(i)
            transcriptIDs.append(transcript)
            transcriptType.append(metric_transcript_bioType(transcript))
            i += 1
    else:
        transcriptNum.append(i)
        transcriptIDs.append(ensembl_transcriptIDs)
        transcriptType.append(metric_transcript_bioType(ensembl_transcriptIDs))
    
    transcriptMetrics = {"TranscriptNum":transcriptNum, "TranscriptID":transcriptIDs, "Type":transcriptType}


    # Generating pd.Dataframe
    return pd.DataFrame(transcriptMetrics)


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
    initialTime = time.perf_counter()
    from paths import *
    main()
    print(time.perf_counter() - initialTime)
if __name__ != "__main__":
    sys.path.append('/home/user1/Dropbox/Research/Neurobiology_PhD/Rotations/Huang/Projects/CellReadR/Code') 
    from kCellReadR.paths import *


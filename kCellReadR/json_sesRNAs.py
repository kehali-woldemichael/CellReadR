import pandas as pd 
import json 
import sys
import time

def json_output_sesRNAs():
    '''Ouputs pandas.DataFrame as json ... for node.js Express server
    sesRNAs information (bioType)'''

    # Reading in necessary parameters 
    # Gene parameters
    speciesName = sys.argv[1]
    geneName = sys.argv[2]
    spliceVariant = int(sys.argv[3])
    searchSeq = sys.argv[4]
    variantTable = sys.argv[5]
    # sesRNA parameters
    seqDirection = sys.argv[6]
    len_sesRNA = int(sys.argv[7])
    minTGG = int(sys.argv[8])
    maxStop = int(sys.argv[9])
    choiceATG = sys.argv[10]

    minGC = int(sys.argv[11])
    maxGC = int(sys.argv[12])
    dist_cTGG = int(sys.argv[13])
    dist_stop_cTGG = int(sys.argv[14])

    ensembl_transcriptIDs = return_ensemblTranscriptIDs(speciesName, geneName)
    variantTable = table_transcriptsInfo(ensembl_transcriptIDs)

    rC_exon_records, C_exon_records, CDS, cDNA = load_referenceSequences(speciesName, geneName, spliceVariant)
    parameters = parameters_sesRNA(speciesName, geneName,  spliceVariant, seqDirection, len_sesRNA, minTGG, maxStop, choiceATG, minGC, maxGC, dist_cTGG, dist_stop_cTGG)
    # parameters.print_parameters()

    if searchSeq == 'CDS': chosen_searchSeq = CDS
    elif searchSeq == 'cDNA': chosen_searchSeq = cDNA

    all_sesRNAs, all_sequenceMetrics, all_sesRNA_objs = generate_all_sesRNAs(rC_exon_records, C_exon_records, chosen_searchSeq, parameters, variantTable)

    # Generating pd.Dataframe
    df = pd.DataFrame(all_sequenceMetrics)
    # Converting DataFrame to json and dumping it to std.out
    df_json = df.reset_index().to_json(orient="values")
    parsed = jsons.loads(df_json)
    print(jsons.dumps(parsed, indent=4))
    sys.stdout.flush()

if __name__ == "__main__":
    from ensembl import *
    from sesRNAs import *
    from IO import *
    from paths import *
    json_output_sesRNAs(_)

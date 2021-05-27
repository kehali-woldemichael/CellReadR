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
    spliceVariant = sys.argv[3]
    # sesRNA parameters
    seqDirection = sys.argv[4]
    len_sesRNA = sys.argv[5]
    minTGG = sys.argv[6]
    maxStop = sys.argv[7]
    choiceATG = sys.argv[8]

    minGC = sys.argv[9]
    maxGC = sys.argv[10]
    dist_cTGG = sys.argv[11]
    dist_stop_cTGG = sys.argv[12]

    # Loading necessary gene sequences 
    rC_exon_records, C_exon_records, CDS, cDNA = load_referenceSequences(geneName, speciesName, spliceVariant)
    # Defining parameters 
    parameters = parameters_sesRNA(speciesName, geneName, 
                                    seqDirection, isoForm, 
                                    len_sesRNA, minTGG, maxStop, minGC, maxGC, 
                                    dist_cTGG, dist_stop_cTGG, choiceATG)

    # Generating pd.Dataframe
    df = pd.DataFrame(transcriptMetrics)
    # Converting DataFrame to json and dumping it to std.out
    df_json = df.reset_index().to_json(orient="values")
    parsed = jsons.loads(df_json)
    print(jsons.dumps(parsed, indent=4))

if __name__ == "__main__":
    from ensembl import *
    from paths import *
    json_output_sesRNAs()

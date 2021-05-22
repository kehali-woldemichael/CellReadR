import pandas as pd 
import json 
import sys
import time

def json_output_spliceVariants():
    '''Ouputs pandas.DataFrame as json ... for node.js Express server
    Splice variant information (bioType)'''

    speciesName = sys.argv[1]
    geneName = sys.argv[2]

    # Generating transcriptIDs for gene
    ensembl_transcriptIDs = return_ensemblTranscriptIDs(speciesName, geneName)
    # Generating transcript splice variant metrics 
    transcriptMetrics = table_transcripsInfo(ensembl_transcriptIDs)
    # Generating pd.Dataframe
    df = pd.DataFrame(transcriptMetrics)

    # Converting DataFrame to json and dumping it to std.out
    df_json = df.reset_index().to_json(orient="values")
    parsed = jsons.loads(df_json)
    print(jsons.dumps(parsed, indent=4))

if __name__ == "__main__":
    from ensembl import *
    from paths import *
    json_output_spliceVariants()

import pandas as pd
import numpy as np
import jsons

# Generating randow dataframe for testing of transfer and rending with node/react
df = pd.DataFrame(np.random.randint(0,100,size=(15, 4)), columns=list('ABCD'))
# Uses indecies as IDs for first row 
df.index.name = 'ID'
df.reset_index(inplace=True)

# Converting DataFrame to json and dumping it to std.out
df_json = df.to_json(orient="split")
parsed = jsons.loads(df_json)
print(jsons.dumps(parsed, indent=4))


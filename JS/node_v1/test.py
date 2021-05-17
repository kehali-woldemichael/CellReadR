import pandas as pd
import numpy as np

df = pd.DataFrame(np.random.randint(0,100,size=(15, 4)), columns=list('ABCD'))
print(df.head())

df_json = df.to_json(orient="split")
parsed = json.loads(df_json)
json.dumps(parsed, indent=4)


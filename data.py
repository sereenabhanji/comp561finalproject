import pandas as pd
import numpy as np

df = pd.read_csv("/Users/sereenabhanji/Downloads/clinvar_result.txt", sep="\t")
arr = df.to_numpy()
num_vars = len(arr)
drop_rows = []

# col 4 describes pathogenicity
for i in range(num_vars):
    curr = arr[i][4]
    if "Pathogen" not in arr[i][4]:
        drop_rows.append(i)

df = df.drop(labels=drop_rows, axis=0)
arr = df.to_numpy()

print("hello")
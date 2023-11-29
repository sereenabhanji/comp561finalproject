import pandas as pd
import numpy as np
from requests import get, codes as http_code
import ratelimit
from typing import Any


def query_by_rsid(rs_id):
    # this function is adapted from the ncbi's dbsnp
    access = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + str(rs_id[2:]) + "/frequency"
    reply = get(access)
    if reply.status_code != http_code.ok:
        return
        # raise Exception("Request failed: {}\n{}".format(
            # reply.status_code, access))

    content_type = reply.headers['content-type']
    if content_type != 'application/json':
        raise Exception("Unexpected content type: {}\n{}".format(
            content_type, access))

    return reply


df = pd.read_csv("/Users/sereenabhanji/Downloads/clinvar_result.txt", sep="\t")
arr = df.to_numpy()
num_vars = len(arr)
drop_rows = set()

# col 4 describes pathogenicity
for i in range(num_vars):
    curr = arr[i][4]
    if "Pathogen" not in arr[i][4]:
        drop_rows.add(i)
    elif arr[i][13] is np.nan:
        drop_rows.add(i)
    elif "Cystic fibrosis" not in arr[i][3]:
        drop_rows.add(i)
    elif "no assertion criteria provided" in arr[i][5]:
        drop_rows.add(i)
    if type(arr[i][13]) is str:
        if '|' in arr[i][13]:
            rs = ""
            j = 0
            s = arr[i][13][j]
            while s != "|":
                rs += s
                j += 1
                s = arr[i][13][j]

            arr[i][13] = rs


drop_rows = list(drop_rows)
df = df.drop(labels=drop_rows, axis=0)
arr = df.to_numpy()

'''
alfa_data = []
for snp in arr:
    alfa_data.append(query_by_rsid(snp[13])) '''


print("hello")
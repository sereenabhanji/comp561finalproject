# this code is adapted from the NLM jupyter notebook
# found at https://github.com/ncbi/dbsnp/blob/master/tutorials/Variation%20Services/Jupyter_Notebook/by_rsid.ipynb

import requests as re
import ratelimit as rl

reply = re.get("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{}/frequency".format(16))
reply.json()


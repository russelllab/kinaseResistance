#data_reformatter

import pickle
import sys

mutations = sys.argv[1]
features  = sys.argv[2]
with open(mutations,"rb") as f:
    data_align = pickle.load(f)

with open("Mutational_Infofile.txt", 'w') as out:
    print(data_align, file= out)

with open(features,"rb") as ff:
    data_feat = pickle.load(ff)

with open("Features_Infofile.txt", 'w') as out:
    print(data_feat, file= out)


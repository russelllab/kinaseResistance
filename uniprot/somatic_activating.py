import sys
import gzip
import csv
import pandas as pd
import numpy as np


if __name__ == '__main__':
    file_path = {
        "uniprot": "../uniprot/uniprot_muts_activark_with_annotations.txt.gz",
        "cosmic": "../cosmic/mutation_counts_all.tsv.gz",
        "mut_mined_checked": "../data/mutations/final_mined_RR_checked_checked-again.txt.gz",
        "mut_mined_all": "../data/mutations/medline_muts_21_Dec_2022_pooled_mp2_cleaned.txt.gz"
    }

    df_uni = pd.read_csv(file_path["uniprot"], sep="\t")
    df_cosmic = pd.read_csv(file_path["cosmic"], sep="\t", header=None, names=["UserInput", "Cosmic"])

    df_rob = []
    with gzip.open(file_path["rob_muts"], "rt") as f:
        for i, row in enumerate(csv.reader(f, delimiter="\t")):
            if i > 0:
                df_rob.append({
                    "UserInput": row[0],
                    "GeneName": row[1],
                    "KnownADR": row[2],
                    "Pubmed": row[3]
                })
    df_rob = pd.DataFrame(df_rob)
    df_rob["n_PMIDs"] = df_rob["Pubmed"].str.count(",") + 1
    df_rob["n_PMIDs"] = df_rob["n_PMIDs"].replace(np.inf, 0).replace(np.nan, 0)
    df_rob["n_PMIDs"] = df_rob["n_PMIDs"].astype('int')

    # df_rob2 = []
    # with gzip.open(file_path["rob_muts_all"], "rt") as f:
    #     for i, row in enumerate(csv.reader(f, delimiter="\t")):
    #         if i > 0:
    #             df_rob2.append({
    #                 "UserInput": row[0],
    #                 "GeneName":  row[2].split("/")[0],
    #                 "n_PMIDS": row[5]
    #             })
    # df_rob2 = pd.DataFrame(df_rob2)
    # print(df_rob2.head())
    # sys.exit()


    df = pd.concat( [df_uni[["UserInput", "GeneName", "KnownADR", "n_PMIDs","Note", "Type"]],
                     df_rob[["UserInput", "GeneName", "KnownADR", "n_PMIDs"]]],
                    ignore_index=True
                    ).drop_duplicates(["UserInput"], keep="last")

    df = df.set_index("UserInput").join(df_cosmic.set_index("UserInput")["Cosmic"], on="UserInput", how="left")
    df["Cosmic"] = df["Cosmic"].fillna(0).astype('int')


    df.to_csv("merged_uniprot_rob2.tsv", sep="\t")

    print( df[(df["KnownADR"].str.contains("activ") | df["KnownADR"].str.contains("increase"))
              & (df["Cosmic"]>0)
              & (df["n_PMIDs"]>0)].sort_values("Cosmic")[["GeneName", "KnownADR", "n_PMIDs", "Cosmic", "Type", "Note"]] )


    sys.exit()
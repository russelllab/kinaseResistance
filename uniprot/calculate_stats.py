import sys
import os
import gzip
import csv
import pandas as pd
import numpy as np


def merge_predictions_with_uniprot_annotations(path):
    """
    Combine Activark predictions with UniProt annotations.
    """
    if os.path.exists(path["merged_file"]):
        return pd.read_csv(path["merged_file"], sep="\t", compression="gzip")
    else:
        df_act = pd.read_csv(path["activark"], sep=r"\s+", comment='#')
        df_uni = pd.read_csv(path["uniprot"], sep="\t", comment='#')
        df_uni.rename(columns={"UniProt/Mutation": "UserInput"}, inplace=True)
        df = pd.merge(df_act, df_uni[["UserInput", "Type", "Note", "Pubmed"]],
                      on="UserInput", how="left")
        df["n_PMIDs"] = df["Pubmed"].str.count(",") + 1
        df["n_PMIDs"] = df["n_PMIDs"].replace(np.inf, 0).replace(np.nan, 0)
        df["n_PMIDs"] = df["n_PMIDs"].astype('int')
        df.to_csv(path["merged_file"], sep="\t", index=False, compression="gzip")
    return df


def compute_stats(df, output_file):
    """
    Compute statistics for the given dataframe.
    """
    final_df = [{
        "Set": "ALL",
        "Total_Vars": str(len(df.index)),
        "Total_Prots": str(len(df["UniProtAcc"].unique())),

        "KnownA_Vars": str(
            len(df[df["KnownADR"].str.contains("activating") | df["KnownADR"].str.contains("increase")].index)),
        "KnownA_%": f"{round(len(df[df['KnownADR'].str.contains('activating') | df['KnownADR'].str.contains('increase')].index) / len(df.index) * 100, 2)}%",
        "KnownR_Vars": str(len(df[df["KnownADR"].str.contains("resistance")].index)),
        "KnownR_%": f"{round(len(df[df['KnownADR'].str.contains('resistance')].index) / len(df.index) * 100, 2)}%",
        "KnownD_Vars": str(
            len(df[df["KnownADR"].str.contains("decrease") | df["KnownADR"].str.contains("loss")].index)),
        "KnownD_%": f"{round(len(df[df['KnownADR'].str.contains('decrease') | df['KnownADR'].str.contains('loss')].index) / len(df.index) * 100, 2)}%",
        "KnownN_Vars": str(len(df[df["KnownADR"].str.contains("neutral")].index)),
        "KnownN_%": f"{round(len(df[df['KnownADR'].str.contains('neutral')].index) / len(df.index) * 100, 2)}%",

        "PredA_Vars": str(len(df[df["AIvLD"] > 0.5].index)),
        "PredA_%": f"{round(len(df[df['AIvLD'] > 0.5].index) / len(df.index) * 100, 2)}%",
        "PredD_Vars": str(len(df[df["AIvLD"] <= 0.5].index)),
        "PredD_%": f"{round(len(df[df['AIvLD'] <= 0.5].index) / len(df.index) * 100, 2)}%",
        "PredR_Vars": str(len(df[df["RvN"] > 0.5].index)),
        "PredR_%": f"{round(len(df[df['RvN'] > 0.5].index) / len(df.index) * 100, 2)}%",
    }]

    for type in ["MUTAGEN", "VARIANT"]:
        df2 = df[df["Type"] == type]
        d ={
            "Set": type,
            "Total_Vars": str(len(df2.index)),
            "Total_Prots": str(len(df2["UniProtAcc"].unique())),

            "KnownA_Vars": str(len(df2[df2["KnownADR"].str.contains("activating") | df2["KnownADR"].str.contains("increase")].index)),
            "KnownA_%": f"{round(len(df2[df2['KnownADR'].str.contains('activating') | df2['KnownADR'].str.contains('increase')].index) / len(df2.index) * 100, 2)}%",
            "KnownR_Vars": str(len(df2[df2["KnownADR"].str.contains("resistance")].index)),
            "KnownR_%": f"{round(len(df2[df2['KnownADR'].str.contains('resistance')].index) / len(df2.index) * 100, 2)}%",
            "KnownD_Vars": str(len(df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")].index)),
            "KnownD_%": f"{round(len(df2[df2['KnownADR'].str.contains('decrease') | df2['KnownADR'].str.contains('loss')].index) / len(df2.index) * 100, 2)}%",
            "KnownN_Vars": str(len(df2[df2["KnownADR"].str.contains("neutral")].index)),
            "KnownN_%": f"{round(len(df2[df2['KnownADR'].str.contains('neutral')].index) / len(df2.index) * 100, 2)}%",
        }
        df_a = df2[(df2["A"] > df2["D"]) & (df2["A"] > df2["N"]) & (df2["AIvLD"] > 0.5)]
        df_d = df2[(df2["D"] > df2["A"]) & (df2["D"] > df2["N"]) % (df2["AIvLD"] <= 0.5)]
        df_n = df2[(df2["N"] > df2["A"]) & (df2["N"] > df2["D"])]
        df_r = df2[df2["RvN"] > 0.5]

        d["PredA_Vars"] = str(len(df_a.index))
        d["PredA_%"] = f"{round(len(df_a.index) / len(df2.index) * 100, 2)}%"
        d["PredD_Vars"] = str(len(df_d.index))
        d["PredD_%"] = f"{round(len(df_d.index) / len(df2.index) * 100, 2)}%"
        d["PredN_Vars"] = str(len(df_n.index))
        d["PredN_%"] = f"{round(len(df_n.index) / len(df2.index) * 100, 2)}%"
        d["PredNA_Vars"] = str(len(df_n[df_n["AIvLD"] > 0.5].index))
        d["PredNA_%"] = f"{round(len(df_n[df_n['AIvLD'] > 0.5].index) / len(df2.index) * 100, 2)}%"
        d["PredND_Vars"] = str(len(df_n[df_n["AIvLD"] <= 0.5].index))
        d["PredND_%"] = f"{round(len(df_n[df_n['AIvLD'] <= 0.5].index) / len(df2.index) * 100, 2)}%"
        d["PredR_Vars"] = str(len(df_r.index))
        d["PredR_%"] = f"{round(len(df_r.index) / len(df2.index) * 100, 2)}%"

        final_df.append(d)

    final_df = pd.DataFrame(final_df)
    final_df.to_csv(output_file, sep="\t", index=False)


if __name__ == '__main__':
    file_path = {
        "activark": "uniprot_muts_activark.txt.gz",
        "uniprot": "allUniProtVariants.txt",
        "merged_file": "uniprot_muts_activark_with_annotations.txt.gz",
        "stats": "stats_from_uniprot_activark.txt",
    }

    df = merge_predictions_with_uniprot_annotations(file_path)

    compute_stats(df, file_path["stats"])

import gzip
import pandas as pd


if __name__ == '__main__':
    file_path = {
        "activark": "uniprot_muts_activark.txt.gz",
        "uniprot": "allUniProtVariants.txt",
        "merged_file": "activark_uniprot_merged.tsv",
        "stats": "stats_from_uniprot_activark.txt"
    }

    df_act = pd.read_csv(file_path["activark"], sep=r"\s+", comment='#')
    df_uni = pd.read_csv(file_path["uniprot"], sep="\t", comment='#')
    df_uni.rename(columns={"UniProt/Mutation": "UserInput"}, inplace=True)
    df = pd.merge(df_act, df_uni[["UserInput", "Type", "Note", "Pubmed"]], 
                  on="UserInput", how="left")
    df.to_csv(file_path["merged_file"], sep="\t", index=False)

    final_df = []
    final_df.append({
        "Set": "ALL",
        "Total_Vars": str(len(df.index)),
        "Total_Prots": str(len(df["UniProtAcc"].unique())),

        "KnownA_Vars": str(len(df[df["KnownADR"].str.contains("activating") | df["KnownADR"].str.contains("increase")].index)),
        "KnownA_%": f"{round(len(df[df['KnownADR'].str.contains('activating') | df['KnownADR'].str.contains('increase')].index) / len(df.index) * 100, 2)}%",
        "KnownR_Vars": str(len(df[df["KnownADR"].str.contains("resistance")].index)),
        "KnownR_%": f"{round(len(df[df['KnownADR'].str.contains('resistance')].index) / len(df.index) * 100, 2)}%",
        "KnownD_Vars": str(len(df[df["KnownADR"].str.contains("decrease") | df["KnownADR"].str.contains("loss")].index)),
        "KnownD_%": f"{round(len(df[df['KnownADR'].str.contains('decrease') | df['KnownADR'].str.contains('loss')].index) / len(df.index) * 100, 2)}%",
        "KnownN_Vars": str(len(df[df["KnownADR"].str.contains("neutral")].index)),
        "KnownN_%": f"{round(len(df[df['KnownADR'].str.contains('neutral')].index) / len(df.index) * 100, 2)}%",

        "PredA_Vars": str(len(df[df["AIvLD"] > 0.5].index)),
        "PredA_%": f"{round(len(df[df['AIvLD'] > 0.5].index) / len(df.index) * 100, 2)}%",
        "PredD_Vars": str(len(df[df["AIvLD"] <= 0.5].index)),
        "PredD_%": f"{round(len(df[df['AIvLD'] <= 0.5].index) / len(df.index) * 100, 2)}%",
        "PredR_Vars": str(len(df[df["RvN"] > 0.5].index)),
        "PredR_%": f"{round(len(df[df['RvN'] > 0.5].index) / len(df.index) * 100, 2)}%",
    })

    for type in ["MUTAGEN", "VARIANT"]:
        df2 = df[df["Type"] == type]
        final_df.append({
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

            "PredA_Vars": str(len(df2[df2["AIvLD"] > 0.5].index)),
            "PredA_%": f"{round(len(df2[df2['AIvLD'] > 0.5].index) / len(df2.index) * 100, 2)}%",
            "PredD_Vars": str(len(df2[df2["AIvLD"] <= 0.5].index)),
            "PredD_%": f"{round(len(df2[df2['AIvLD'] <= 0.5].index) / len(df2.index) * 100, 2)}%",
            "PredR_Vars": str(len(df2[df2["RvN"] > 0.5].index)),
            "PredR_%": f"{round(len(df2[df2['RvN'] > 0.5].index) / len(df2.index) * 100, 2)}%",
        })

    final_df = pd.DataFrame(final_df)
    final_df.to_csv(file_path["stats"], sep="\t", index=False)
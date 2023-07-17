import sys
import os
import gzip
import csv
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def merge_predictions_with_uniprot_annotations(path):
    """
    Combine Activark predictions with UniProt annotations.
    """
    # if os.path.exists(path["merged_file"]):
    #     return pd.read_csv(path["merged_file"], sep="\t", compression="gzip")
    # else:
    df_act = pd.read_csv(path["activark"], sep=r"\s+", comment='#')
    df_uni = pd.read_csv(path["uniprot"], sep="\t", comment='#')
    df_uni.rename(columns={"UniProt/Mutation": "UserInput"}, inplace=True)
    df = pd.merge(df_act, df_uni[["UserInput", "Type", "Note", "Pubmed"]],
                    on="UserInput", how="left")
    df["n_PMIDs"] = df["Pubmed"].str.count(",") + 1
    df["n_PMIDs"] = df["n_PMIDs"].replace(np.inf, 0).replace(np.nan, 0)
    df["n_PMIDs"] = df["n_PMIDs"].astype('int')
    df["Note"] = df["Note"].str.replace("T(-)B(+)NK(-) SCID", "TBNKSCID", regex=False)
    df["Disease"] = df["Note"].str.extract(r'[^\w]in ([A-Z0-9]{2,})')
    df["Disease"] = df["Disease"].str.replace("TBNKSCID", "T(-)B(+)NK(-) SCID", regex=False)
    df["Note"] = df["Note"].str.replace("TBNKSCID", "T(-)B(+)NK(-) SCID", regex=False)
    for val in ["BC", "CRC", "G1", "GASC", "GIST", "GISTPS", "GLM", "HCC", 
                "LNCR", "MAPK15", "MAP2K1", "MTC", "MEN2A", "MEN2B", "NHL", "OC", 
                "RCCP", "RNA", "SH2B2", "STAP2", "STAT3", "TGCT"]:
        df["Disease"] = df["Disease"].replace(val, np.nan)

    df.sort_values(by=['Disease']).to_csv(path["merged_file"], sep="\t", index=False, compression="gzip")
    # with open("see.tsv", "wt") as f:
    #     f.write("\n".join([str(dis) for dis in sorted(list(set([str(x) for x in df["Disease"].tolist()])))]))
    df[df["Note"].str.contains("in ")]["Note"].to_csv("see.tsv", sep="\t", index=False)
    return df


def compute_stats(df, output_file):
    """
    Compute statistics for the given dataframe.
    """
    final_df = []

    for type in ["ALL", "MUTAGEN", "VARIANT", "DISEASE"]:
        if type == "ALL":
            df2 = df[df["N"].notnull()]
        elif type == "DISEASE":
            df2 = df[(df["N"].notnull()) & (df["Disease"].notnull())]
        else:
            df2 = df[(df["N"].notnull()) & (df["Type"] == type)]

        known_act = df2[df2["KnownADR"].str.contains("activating") | df2["KnownADR"].str.contains("increase")]
        known_deac = df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")]
        known_neu = df2[df2["KnownADR"].str.contains("neutral")]
        known_res = df2[df2["KnownADR"].str.contains("resistance")]

        df2 = df2[df2["KnownADR"] == "-"]
        df_a = df2[(df2["A"] > df2["D"]) & (df2["A"] > df2["N"]) & (df2["AIvLD"] > 0.5)]
        df_d = df2[(df2["D"] > df2["A"]) & (df2["D"] > df2["N"]) & (df2["AIvLD"] <= 0.5)]
        df_n = df2[(df2["N"] > df2["A"]) & (df2["N"] > df2["D"])]
        df_r = df2[df2["RvN"] > 0.5]

        d ={
            "Set": type,
            "Total_Vars": str(len(df2.index)),
            "Total_Prots": str(len(df2["UniProtAcc"].unique())),

            "KnownA_Vars": str(len(known_act.index)),
            "KnownA_%": f"{round(len(known_act.index) / len(df2.index) * 100, 2)}%",
            "KnownD_Vars": str(len(known_deac.index)),
            "KnownD_%": f"{round(len(known_deac.index) / len(df2.index) * 100, 2)}%",
            "KnownN_Vars": str(len(known_neu.index)),
            "KnownN_%": f"{round(len(known_neu.index) / len(df2.index) * 100, 2)}%",
            "KnownR_Vars": str(len(known_res.index)),
            "KnownR_%": f"{round(len(known_res.index) / len(df2.index) * 100, 2)}%",

            "PredA_Vars": str(len(df_a.index)),
            "PredA_%": f"{round(len(df_a.index) / len(df2.index) * 100, 2)}%",
            "PredD_Vars": str(len(df_d.index)),
            "PredD_%": f"{round(len(df_d.index) / len(df2.index) * 100, 2)}%",
            "PredN_Vars": str(len(df_n.index)),
            "PredN_%": f"{round(len(df_n.index) / len(df2.index) * 100, 2)}%",
            "PredNA_Vars": str(len(df2[df2["AIvLD"] > 0.5].index)),
            "PredNA_%": f"{round(len(df2[df2['AIvLD'] > 0.5].index) / len(df2.index) * 100, 2)}%",
            "PredND_Vars": str(len(df2[df2["AIvLD"] <= 0.5].index)),
            "PredND_%": f"{round(len(df2[df2['AIvLD'] <= 0.5].index) / len(df2.index) * 100, 2)}%",
            "PredR_Vars": str(len(df_r.index)),
            "PredR_%": f"{round(len(df_r.index) / len(df2.index) * 100, 2)}%"
        }

        final_df.append(d)

    final_df = pd.DataFrame(final_df)
    final_df.to_csv(output_file, sep="\t", index=False)

    return


def compute_stats_for_diseases(df, output_file):
       
    df["Pred"] = np.where(df["AIvLD"] > 0.5, "A", "D")
    # df["Pred"] = df.apply(lambda row: decide_prediction(row), axis=1)
    # df = df[df["Pred"] != "-"]
    df = df[(df["N"].notnull()) & (df["Disease"].notnull()) & (df["Type"] == "VARIANT")]
    #df = df[df["Disease"].str.contains("NS")] # Find a particular disease

    final_df = []
    for disease in df["Disease"].unique():
        df3 = df[df["Disease"] == disease]
        unk = df3[df3["KnownADR"] == "-"]
        d = {
            "Disease": disease, 
            "VarNumber": len(df3.index),
            "GeneNumber": len(df3["GeneName"].unique()),
            "GeneList": ", ".join(list(df3["GeneName"].unique())),
            "KnownA": len(df3[df3["KnownADR"].str.contains("activating") | df3["KnownADR"].str.contains("increase")].index),
            "KnownD": len(df3[df3["KnownADR"].str.contains("decrease") | df3["KnownADR"].str.contains("loss")].index),
            "KnownR": len(df3[df3["KnownADR"].str.contains("resistance")].index),

            "PredA": len(df3[df3["Pred"] == "A"].index),
            "PredD": len(df3[df3["Pred"] == "D"].index),
            "PredA%": round(len(df3[df3['Pred'] == 'A'].index) / len(df3.index), 2),
            "PredD%": round(len(df3[df3['Pred'] == 'D'].index) / len(df3.index), 2),

            "UnkPredA": len(unk[unk["Pred"] == "A"].index),
            "UnkPredD": len(unk[unk["Pred"] == "D"].index),
        }
        try:
            d["UnkPredA%"] = round(len(unk[unk['Pred'] == 'A'].index) / len(unk.index), 2)
        except:
            d["UnkPredA%"] = 0
        try:
            d["UnkPredD%"] = round(len(unk[unk['Pred'] == 'D'].index) / len(unk.index), 2)
        except:
            d["UnkPredD%"] = 0

        if d["KnownA"] > 0 and d["KnownD"] == 0:
            d["KnownVariants"] = "Only activating"
        elif d["KnownD"] > 0 and d["KnownA"] == 0:
            d["KnownVariants"] = "Only deactivating"
        elif d["KnownA"] == 0 and d["KnownD"] == 0:
            d["KnownVariants"] = "None"
        else:
            d["KnownVariants"] = "Both"

        # try:
        #     d["known"] = d["KnownA"]/(d["KnownA"]+d["KnownD"]) - d["KnownD"]/(d["KnownA"]+d["KnownD"])
        # except:
        #     d["known"] = 0

        final_df.append(d)

    df_all = pd.DataFrame(final_df)
    df_all.sort_values(by=['Disease']).to_csv(output_file.replace(".tsv", "_all.tsv"), sep="\t", index=False)
    
    df_unk = df_all[df_all["VarNumber"] > (df_all["KnownA"]+df_all["KnownD"])]
    df_unk.sort_values(by=['Disease']).to_csv(output_file.replace(".tsv", "_unkvars.tsv"), sep="\t", index=False)

    sns.set_theme(style="whitegrid")
    col_order = ["Only activating", "Only deactivating", "Both", "None"]
    palette = {"Only activating": "#28b463", "Only deactivating": "#cb4335", "None": "#839192", "Both": "#f1c40f"}

    ## For all diseases  
    g = sns.relplot(data=df_all, x="PredA%", y="PredD%", 
                    size="VarNumber", hue="KnownVariants", sizes=(30, 300), alpha=1, height=6, aspect=1.3,
                    palette=palette)
    g.despine(left=True, bottom=True)
    g.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (all)')
    label_point(plt.gca(), df_all["PredA%"], df_all["PredD%"], df_all["Disease"], df_all["KnownVariants"],
                col_order) 
    g.savefig("disease_plot_all.pdf", dpi=300, format="pdf")

    g1 = sns.relplot(data=df_all, x="PredA%", y="PredD%", col="KnownVariants", col_order=col_order,
                     size="VarNumber", hue="KnownVariants", sizes=(30, 300), alpha=1, height=5, aspect=1,
                     palette=palette)
    g1.despine(left=True, bottom=True)
    g1.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (all)')
    label_point_split(plt.gcf().get_axes(), df_all["PredA%"], df_all["PredD%"], df_all["Disease"], df_all["KnownVariants"],
                col_order) 
    g1.savefig("disease_plot_all_split.pdf", dpi=300, format="pdf")

    # For diseases with unKnownVariants
    g2 = sns.relplot(data=df_unk, x="UnkPredA%", y="UnkPredD%", 
                     size="VarNumber", hue="KnownVariants", sizes=(30, 300), alpha=1, height=6, aspect=1.3,
                     palette=palette)
    g2.despine(left=True, bottom=True)
    g2.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (all)')
    label_point(plt.gca(), df_unk["UnkPredA%"], df_unk["UnkPredD%"], df_unk["Disease"], df_unk["KnownVariants"],
                col_order) 
    g2.savefig("disease_plot_unkvars.pdf", dpi=300, format="pdf")

    g3 = sns.relplot(data=df_unk, x="UnkPredA%", y="UnkPredD%", col="KnownVariants", col_order=col_order,
                    size="VarNumber", hue="KnownVariants", sizes=(30, 300), alpha=1, height=5, aspect=1,
                    palette=palette)
    g3.despine(left=True, bottom=True)
    g3.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (with unk)')
    label_point_split(plt.gcf().get_axes(), df_unk["UnkPredA%"], df_unk["UnkPredD%"], df_unk["Disease"], df_unk["KnownVariants"],
                col_order)  
    g3.savefig("disease_plot_unkvars_split.pdf", dpi=300, format="pdf")

    return


def decide_prediction(row):
    if row["A"] > 0.5:
        return "A"
    elif row["D"] > 0.5:
        return "D"
    else:
        return "-"


def label_point(ax, x, y, val, type, col_order):
    color = {"Only activating": "#28b463", "Only deactivating": "#cb4335", "None": "#839192", "Both": "#f1c40f"}
    a = pd.concat({'x': x, 'y': y, 'label': val, "type":type}, axis=1)
    for cat in col_order:
        for i, point in a[a["type"]==cat].iterrows():
            if point["type"] == "Only deactivating":
                ax.text(point['x']-.10, point['y']-.02, str(point['label']), color=color[point["type"]])
            else:
                ax.text(point['x']+.02, point['y'], str(point['label']), color=color[point["type"]])


def label_point_split(axes, x, y, val, type, col_order):
    color = {"Only activating": "#28b463", "Only deactivating": "#cb4335", "None": "#839192", "Both": "#f1c40f"}
    a = pd.concat({'x': x, 'y': y, 'label': val, "type":type}, axis=1)
    for ax, cat in zip(axes, col_order):
        for i, point in a[a["type"]==cat].iterrows():
            # if point["type"] == "Only deactivating":
            #     ax.text(point['x']-.10, point['y']-.02, str(point['label']), color=color[point["type"]])
            # else:
            ax.text(point['x']+.02, point['y'], str(point['label']), color=color[point["type"]])


if __name__ == '__main__':
    file_path = {
        "activark": "uniprot_muts_activark.txt.gz",
        "uniprot": "allUniProtVariants.txt",
        "merged_file": "uniprot_muts_activark_with_annotations.txt.gz",
        "stats": "stats_from_uniprot_activark.tsv",
        "stats_dis": "stats_uniprot_diseases.tsv"
    }

    df = merge_predictions_with_uniprot_annotations(file_path)
 
    compute_stats(df, file_path["stats"])

    compute_stats_for_diseases(df, file_path["stats_dis"])

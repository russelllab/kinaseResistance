import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys


def merge_predictions_with_uniprot_annotations(path):
    """
    Combine Activark predictions with UniProt annotations.
    """
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

    return df


def decide_AvD_prediction(row):
    """
    Make decision on whether variant is predicted to be activating/deactivating
    with different confidence levels, or uncertain.
    """
   
    if row["AIvLD"] > 0.7 and row["A"] > row["D"] and row["A"] > row["N"]:
        return "Activating (High)"
    elif row["AIvLD"] > 0.5 and row["A"] > row["D"] and row["A"] > row["N"]:
        return "Activating (Medium)"
    elif row["AIvLD"] > 0.5 and row["N"] > row["A"] > row["D"]:
        return "Activating (Low)"
    elif row["AIvLD"] < 0.3 and row["D"] > row["A"] and row["D"] > row["N"]:
        return "Deactivating (High)"
    elif row["AIvLD"] < 0.5 and row["D"] > row["A"] and row["D"] > row["N"]:
        return "Deactivating (Medium)"
    elif row["AIvLD"] < 0.5 and row["N"] > row["D"] > row["A"]:
        return "Deactivating (Low)"
    else:
        return "Uncertain"
    

def compute_stats_for_variants(df, file_path):
    """
    Compute some numbers for the variants in dataframe.
    """

    df["AD_prediction"] = df.apply(decide_AvD_prediction, axis=1)
    
    final_df = []
    for type in ["ALL", "MUTAGEN", "VARIANT", "DISEASE"]:
        if type == "ALL":
            df2 = df[df["N"].notnull()]
        elif type == "MUTAGEN" or type == "VARIANT":
            df2 = df[(df["N"].notnull()) & (df["Type"] == type)]
        elif type == "DISEASE":
            df2 = df[(df["N"].notnull()) & (df["Type"] == "VARIANT") & (df["Disease"].notnull())]

        known_act = df2[df2["KnownADR"].str.contains("activating") | df2["KnownADR"].str.contains("increase")]
        known_deac = df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")]
        known_neu = df2[df2["KnownADR"].str.contains("neutral")]
        known_res = df2[df2["KnownADR"].str.contains("resistance")]

        pred_act = df2[df2["AD_prediction"].str.contains("Activating")]
        pred_deact = df2[df2["AD_prediction"].str.contains("Deactivating")]
        pred_unce = df2[df2["AD_prediction"].str.contains("Uncertain")]
        pred_res = df2[df2["RvN"] > 0.5]

        d = {
            "VarSet": type,
            "Vars": str(len(df2.index)),
            "Prots": str(len(df2["UniProtAcc"].unique())),

            "KnownA_Vars": str(len(known_act.index)),
            "KnownA_%": f"{round(len(known_act.index) / len(df2.index) * 100, 2)}%",
            "KnownD_Vars": str(len(known_deac.index)),
            "KnownD_%": f"{round(len(known_deac.index) / len(df2.index) * 100, 2)}%",
            "KnownN_Vars": str(len(known_neu.index)),
            "KnownN_%": f"{round(len(known_neu.index) / len(df2.index) * 100, 2)}%",
            "KnownR_Vars": str(len(known_res.index)),
            "KnownR_%": f"{round(len(known_res.index) / len(df2.index) * 100, 2)}%",

            "PredA_Vars": str(len(pred_act.index)),
            "PredA_%": f"{round(len(pred_act.index) / len(df2.index) * 100, 2)}%",
            "PredD_Vars": str(len(pred_deact.index)),
            "PredD_%": f"{round(len(pred_deact.index) / len(df2.index) * 100, 2)}%",
            "PredU_Vars": str(len(pred_unce.index)),
            "PredU_%": f"{round(len(pred_unce.index) / len(df2.index) * 100, 2)}%",
            "PredR_Vars": str(len(pred_res.index)),
            "PredR_%": f"{round(len(pred_res.index) / len(df2.index) * 100, 2)}%",

            "PredR_ActVars": str(len(pred_res[pred_res['AD_prediction'].str.contains('Activating')].index)),
            "PredR_ActVars_%": f"{round(len(pred_res[pred_res['AD_prediction'].str.contains('Activating')].index) / len(pred_res.index) * 100, 2)}%",
            "PredR_DeactVars": str(len(pred_res[pred_res['AD_prediction'].str.contains('Deactivating')].index)),
            "PredR_DeactVars_%": f"{round(len(pred_res[pred_res['AD_prediction'].str.contains('Deactivating')].index) / len(pred_res.index) * 100, 2)}%",
        }

        final_df.append(d)

    final_df = pd.DataFrame(final_df)
    final_df.to_csv(file_path["stats"], sep="\t", index=False)

    columns = {
    "VarSet": "UniProt variant data set: all variants ('ALL'), only variants from mutagenesis experiments ('MUTAGEN')"+
    ", only real variants ('VARIANTS') or only variants associated to a hereditary genetic disease ('DISEASE')",  
    "Vars": "Number of unique variants",
    "Prots": "Number of proteins",
    "KnownA_Vars": "Number of unique known activating variants",
    "KnownA_%": "Fraction from total number of unique variants",
    "KnownD_Vars": "Number of unique known deactivating variants",
    "KnownD_%": "Fraction from total number of unique variants",
    "KnownN_Vars": "Number of unique known neutral variants",
    "KnownN_%": "Fraction from total number of unique variants",
    "KnownR_Vars": "Number of unique known resistance variants",
    "KnownR_%": "Fraction from total number of unique variants",
    "PredA_Vars": "Number of unique variants predicted as activating",
    "PredA_%": "Fraction from total number of unique variants",
    "PredD_Vars": "Number of unique variants predicted as deactivating",
    "PredD_%": "Fraction from total number of unique variants",
    "PredU_Vars": "Number of unique variants predicted as uncertain",
    "PredU_%": "Fraction from total number of unique variants",
    "PredR_Vars": "Number of unique variants predicted as resistance",
    "PredR_%": "Fraction from total number of unique variants",
    "PredR_ActVars": "Number of unique variants predicted as resistance and activating",
    "PredR_ActVars_%": "Fraction from total number of unique variants predicted as resistance",
    "PredR_DeactVars": "Number of unique variants predicted as resistance and deactivating",
    "PredR_DeactVars_%": "Fraction from total number of unique variants predicted as resistance"
    }
    with open(file_path["stats_cols"], "w") as f:
        for k, v in columns.items():
            f.write(f"{k}\t{v}\n")

    return


def compute_stats_for_diseases(df, file_path):
    """
    Compute some numbers for variants in each hereditary disease.
    """

    df["AD_prediction"] = df.apply(decide_AvD_prediction, axis=1)

    df = df[(df["N"].notnull()) & (df["Disease"].notnull()) & (df["Type"] == "VARIANT")]
    # To find a particular disease use:
    # df = df[df["Disease"].str.contains("NS")]
   
    final_df = []
    for disease in df["Disease"].unique():
        df3 = df[df["Disease"] == disease]
        known_act = df3[df3["KnownADR"].str.contains("activating") | df3["KnownADR"].str.contains("increase")]
        known_deac = df3[df3["KnownADR"].str.contains("decrease") | df3["KnownADR"].str.contains("loss")]
        known_neu = df3[df3["KnownADR"].str.contains("neutral")]
        known_res = df3[df3["KnownADR"].str.contains("resistance")]

        pred_act = df3[df3["AD_prediction"].str.contains("Activating")]
        pred_deact = df3[df3["AD_prediction"].str.contains("Deactivating")]
        pred_unce = df3[df3["AD_prediction"].str.contains("Uncertain")]
        pred_res = df3[df3["RvN"] > 0.5]

        unk = df3[df3["KnownADR"] == "-"]


        d = {
            "Disease": disease, 
            "Vars": len(df3.index),
            "Prots": len(df3["GeneName"].unique()),
            "GeneList": ", ".join(list(df3["GeneName"].unique())),
            "KnownA": len(known_act.index),
            "KnownD": len(known_deac.index),
            "KnownR": len(known_res.index),
            "PredA": len(pred_act.index),
            "PredD": len(pred_deact.index),
            "PredU": len(pred_unce.index),
            "UnkPredA": len(unk[unk["AD_prediction"].str.contains("Activating")].index),
            "UnkPredD": len(unk[unk["AD_prediction"].str.contains("Deactivating")].index)
        }
        if d["PredA"]+d["PredD"] == 0:
            continue
        else:
            d["PredA%"] = round(d["PredA"] / (d["PredA"]+d["PredD"]), 2)
            d["PredD%"] = round(d["PredD"] / (d["PredA"]+d["PredD"]), 2)
        

        # Assign category to disease according to known functional variants
        if d["KnownA"]+d["KnownD"] > 1:
            if d["KnownA"] > 0 and d["KnownD"] == 0:
                d["KnownVariants"] = "Only activating"
            elif d["KnownD"] > 0 and d["KnownA"] == 0:
                d["KnownVariants"] = "Only deactivating"
            elif d["KnownA"] == 0 and d["KnownD"] == 0:
                d["KnownVariants"] = "None"
            else:
                d["KnownVariants"] = "Both"
        else:
            d["KnownVariants"] = "None"

        final_df.append(d)


    # Print all diseases with at least 3 variants
    df_all = pd.DataFrame(final_df)
    df_all["UnkPredA%"] = df_all["UnkPredA"] / (df_all["UnkPredA"]+df_all["UnkPredD"])
    df_all["UnkPredD%"] = df_all["UnkPredD"] / (df_all["UnkPredA"]+df_all["UnkPredD"])
    df_all = df_all[df_all["Vars"] > 2]
    df_all.sort_values(by=['Disease']).to_csv(file_path["stats_dis"].replace(".tsv", "_all.tsv"), sep="\t", index=False)
    
    # Print all diseases with at least 1 uncharacterized variant
    df_unk = df_all[df_all["Vars"] > (df_all["KnownA"]+df_all["KnownD"])]
    df_unk.sort_values(by=['Disease']).to_csv(file_path["stats_dis"].replace(".tsv", "_unkvars.tsv"), sep="\t", index=False)

    columns = {
        "Disease": "Disease acronym",
        "Vars": "Number of variants",
        "Prots": "Number of proteins (kinases)",
        "GeneList": "List of proteins (Gene names)",
        "KnownA": "Number of known activating variants",
        "KnownD": "Number of known deactivating variants",
        "KnownR": "Number of known resistance variants",
        "PredA": "Number of variants predicted as activating",
        "PredA%": "Fraction of variants predicted as activating from the total predicted either A or D",
        "PredD": "Number of variants predicted as deactivating",
        "PredD%": "Fraction of variants predicted as deactivating from the total predicted either A or D",
        "PredU": "Number of variants predicted as uncertain",
        "UnkPredA": "Number of uncharacterized variants predicted as activating",
        "UnkPredD": "Number of uncharacterized variants predicted as deactivating",
        "UnkPredA%": "Fraction of uncharacterized variants predicted as activating from the total predicted either A or D",
        "UnkPredD%": "Fraction of uncharacterized variants predicted as deactivating from the total predicted either A or D"
    }
    with open(file_path["stats_dis_cols"], "w") as f:
        for k, v in columns.items():
            f.write(f"{k}\t{v}\n")

    ## Make Plots
    sns.set_theme(style="whitegrid")
    col_order = ["Only activating", "Only deactivating", "Both", "None"]
    palette = {"Only activating": "#28b463", "Only deactivating": "#cb4335", "None": "#839192", "Both": "#f1c40f"}

    ## For all diseases  
    g = sns.relplot(data=df_all, x="PredA%", y="PredD%", 
                    size="Vars", hue="KnownVariants", sizes=(30, 300), alpha=1, height=6, aspect=1.3,
                    palette=palette)
    g.despine(left=True, bottom=True)
    g.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (all)')
    label_point(plt.gca(), df_all["PredA%"], df_all["PredD%"], df_all["Disease"], df_all["KnownVariants"],
                col_order) 
    g.savefig("plots/disease_plot_all.pdf", dpi=300, format="pdf")

    g1 = sns.relplot(data=df_all, x="PredA%", y="PredD%", col="KnownVariants", col_order=col_order,
                     size="Vars", hue="KnownVariants", sizes=(30, 300), alpha=1, height=5, aspect=1,
                     palette=palette)
    g1.despine(left=True, bottom=True)
    g1.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (all)')
    label_point_split(plt.gcf().get_axes(), df_all["PredA%"], df_all["PredD%"], df_all["Disease"], df_all["KnownVariants"],
                col_order) 
    g1.savefig("plots/disease_plot_all_split.pdf", dpi=300, format="pdf")

    # For diseases with unKnownVariants
    g2 = sns.relplot(data=df_unk, x="UnkPredA%", y="UnkPredD%", 
                     size="Vars", hue="KnownVariants", sizes=(30, 300), alpha=1, height=6, aspect=1.3,
                     palette=palette)
    g2.despine(left=True, bottom=True)
    g2.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (all)')
    label_point(plt.gca(), df_unk["UnkPredA%"], df_unk["UnkPredD%"], df_unk["Disease"], df_unk["KnownVariants"],
                col_order) 
    g2.savefig("plots/disease_plot_unkvars.pdf", dpi=300, format="pdf")

    g3 = sns.relplot(data=df_unk, x="UnkPredA%", y="UnkPredD%", col="KnownVariants", col_order=col_order,
                    size="Vars", hue="KnownVariants", sizes=(30, 300), alpha=1, height=5, aspect=1,
                    palette=palette)
    g3.despine(left=True, bottom=True)
    g3.set(xlabel ="Fraction activating", ylabel = "Fraction deactivating")#, title ='Ratio of act/deactivating predicted kinase variants in genetic disease (with unk)')
    label_point_split(plt.gcf().get_axes(), df_unk["UnkPredA%"], df_unk["UnkPredD%"], df_unk["Disease"], df_unk["KnownVariants"],
                col_order)  
    g3.savefig("plots/disease_plot_unkvars_split.pdf", dpi=300, format="pdf")

    return


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
        "uniprot": "data/allUniProtVariants.txt",
        "activark": "ML/uniprot_muts_activark.txt.gz",
        "merged_file": "results/uniprot_muts_activark_with_annotations.txt.gz",
        "stats": "results/stats_from_uniprot_activark.tsv",
        "stats_cols": "results/stats_from_uniprot_activark_cols.tsv",
        "stats_dis": "results/stats_uniprot_diseases.tsv",
        "stats_dis_cols": "results/stats_uniprot_diseases_cols.tsv"
    }

    df = merge_predictions_with_uniprot_annotations(file_path)

    compute_stats_for_variants(df, file_path)

    compute_stats_for_diseases(df, file_path)

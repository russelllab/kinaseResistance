import sys
import pandas as pd
import numpy as np


def get_data_paths():
    """
    Get paths to all data files
    """
    file_paths = {
        "cosmic_ml_all": "results/ML_output_cosmic_all.tsv.gz",
        "cosmic_ml_gws": "results/ML_output_cosmic_gws.tsv.gz",
        "stats": "results/stats_from_cosmic_activark.txt",
        "cols": "results/stats_from_cosmic_activark_cols.txt"
    }
    return file_paths

def decide_AvD_prediction(row):

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


def compute_stats_from_cosmic_ml(file_path):
    """
    Calculate a bunch of stats from the cosmic ML output files. Feel free to add more to it.
    Output is a df printed in a tab-delimited file.
    """
    final_df = []
    for mode, file in zip(["all", "gws"], [file_path["cosmic_ml_all"], file_path["cosmic_ml_gws"]]):
        df = pd.read_csv(file, sep="\t")
        df['GeneRole'] = df['GeneRole'].replace(np.nan, "-")
        df["AD_prediction"] = df.apply(decide_AvD_prediction, axis=1)

        for count_min in [0, 1, 2, 5, 10, 20, 30]:
            df2 = df[df['Cosmic'] > count_min]
            df_onco = df2[df2["GeneRole"].str.contains("oncogene")]
            df_tsg = df2[df2["GeneRole"].str.contains("TSG")]
            df_both = df2[df2["GeneRole"].str.contains("oncogene") & df2["GeneRole"].str.contains("TSG")]
            df_unk = df2[df2["GeneRole"].str.contains("-")]
            known_act = df2[df2["KnownADR"].str.contains("activating") | df2["KnownADR"].str.contains("increase")]
            known_res = df2[df2["KnownADR"].str.contains("resistance")]
            known_deact = df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")]
            known_neut = df2[df2["KnownADR"].str.contains("neutral")]
            pred_act = df2[df2["AD_prediction"].str.contains("Activating")]
            pred_deact = df2[df2["AD_prediction"].str.contains("Deactivating")]
            pred_unce = df2[df2["AD_prediction"].str.contains("Uncertain")]
            pred_res = df2[df2["RvN"] > 0.5]
            d = {
                "VarSet": mode,
                "MinCount": f">{count_min}",
                "Vars": str(len(df2.index)),
                "TotVars": str(df2['Cosmic'].sum()),
                "AvgCount": str(round(df2["Cosmic"].mean(), 2)),
                "Genes": str(df2["GeneName"].nunique()),
                "Oncog": str(df_onco["GeneName"].nunique()), 
                "TSG": str(df_tsg["GeneName"].nunique()),
                "Both": str(df_both["GeneName"].nunique()),
                "Unknown": str(df_unk["GeneName"].nunique()),

                "Oncog_Vars": str(len(df_onco.index)),
                "Oncog_%": f"{round(len(df_onco.index) / len(df2.index) * 100, 2)}%",
                "Oncog_AvgCount": str(round(df_onco["Cosmic"].mean(), 2)),
                "TSG_Vars": str(len(df_tsg.index)),
                "TSG_%": f"{round(len(df_tsg.index) / len(df2.index) * 100, 2)}%",
                "TSG_AvgCount": str(round(df_tsg["Cosmic"].mean(), 2)),

                "KnownA_Vars": str(len(known_act.index)),
                "KnownA_%": f"{round(len(known_act.index) / len(df2.index) * 100, 2)}%",
                "KnownA_TotCount": str(round(known_act["Cosmic"].sum(), 2)),
                "KnownA_TotCount_%": f"{round(known_act['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "KnownA_AvgCount": str(round(known_act["Cosmic"].mean(), 2)),
                "KnownA_Genes": str(known_act["GeneName"].nunique()),
               
                "KnownR_Vars": str(len(known_res.index)),
                "KnownR_%": f"{round(len(known_res.index) / len(df2.index) * 100, 2)}%",
                "KnownR_TotCount": str(round(known_res["Cosmic"].sum(), 2)),
                "KnownR_TotCount_%": f"{round(known_res['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "KnownR_AvgCount": str(round(known_res["Cosmic"].mean(), 2)),
                "KnownR_Genes": str(known_res["GeneName"].nunique()),

                "KnownD_Vars": str(len(known_deact.index)),
                "KnownD_%": f"{round(len(known_deact.index) / len(df2.index) * 100, 2)}%",
                "KnownD_TotCount": str(round(known_deact["Cosmic"].sum(), 2)),
                "KnownD_TotCount_%": f"{round(known_deact['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "KnownD_AvgCount": str(round(known_deact["Cosmic"].mean(), 2)),
                "KnownD_Genes": str(known_deact["GeneName"].nunique()),
                
                "KnownN_Vars": str(len(known_neut.index)),
                "KnownN_%": f"{round(len(known_neut.index) / len(df2.index) * 100, 2)}%",
                "KnownN_TotCount": str(round(known_neut["Cosmic"].sum(), 2)),
                "KnownN_TotCount_%": f"{round(known_neut['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "KnownN_AvgCount": str(round(known_neut["Cosmic"].mean(), 2)),
                "KnownN_Genes": str(known_neut["GeneName"].nunique()),

                "PredA_Vars": str(len(pred_act.index)),
                "PredA_%": f"{round(len(pred_act.index) / len(df2.index) * 100, 2)}%",
                "PredA_TotCount": str(round(pred_act["Cosmic"].sum(), 2)),
                "PredA_TotCount_%": f"{round(pred_act['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "PredA_AvgCount": str(round(pred_act["Cosmic"].mean(), 2)),
                "PredA_Genes": str(pred_act["GeneName"].nunique()),

                "PredD_Vars": str(len(pred_deact.index)),
                "PredD_%": f"{round(len(pred_deact.index) / len(df2.index) * 100, 2)}%",
                "PredD_TotCount": str(round(pred_deact["Cosmic"].sum(), 2)),
                "PredD_TotCount_%": f"{round(pred_deact['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "PredD_AvgCount": str(round(pred_deact["Cosmic"].mean(), 2)),
                "PredD_Genes": str(pred_deact["GeneName"].nunique()),

                "PredU_Vars": str(len(pred_unce.index)),
                "PredU_%": f"{round(len(pred_unce.index) / len(df2.index) * 100, 2)}%",
                "PredU_TotCount": str(round(pred_unce["Cosmic"].sum(), 2)),
                "PredU_TotCount_%": f"{round(pred_unce['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "PredU_AvgCount": str(round(pred_unce["Cosmic"].mean(), 2)),
                "PredU_Genes": str(pred_unce["GeneName"].nunique()),

                "PredR_Vars": str(len(pred_res.index)),
                "PredR_%": f"{round(len(pred_res.index) / len(df2.index) * 100, 2)}%",
                "PredR_TotCount": str(round(pred_res["Cosmic"].sum(), 2)),
                "PredR_TotCount_%": f"{round(pred_res['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "PredR_AvgCount": str(round(pred_res["Cosmic"].mean(), 2)),
                "PredR_Genes": str(pred_res["GeneName"].nunique()),

                "Onco_PredA_Vars": str(len(df_onco[df_onco['AD_prediction'].str.contains('Activating')].index)),
                "Onco_PredA_%": f"{round(len(df_onco[df_onco['AD_prediction'].str.contains('Activating')].index) / len(df_onco.index) * 100, 2)}%",
                "Onco_PredA_AvgCount": str(round(df_onco[df_onco['AD_prediction'].str.contains('Activating')]["Cosmic"].mean(), 2)),
                "Onco_PredD_Vars": str(len(df_onco[df_onco['AD_prediction'].str.contains('Deactivating')].index)),
                "Onco_PredD_%": f"{round(len(df_onco[df_onco['AD_prediction'].str.contains('Deactivating')].index) / len(df_onco.index) * 100, 2)}%",
                "Onco_PredD_AvgCount": str(round(df_onco[df_onco['AD_prediction'].str.contains('Deactivating')]["Cosmic"].mean(), 2)),
                "Onco_PredR_Vars": str(len(df_onco[df_onco["RvN"] > 0.5].index)),
                "Onco_PredR_%": f"{round(len(df_onco[df_onco['RvN'] > 0.5].index) / len(df_onco.index) * 100, 2)}%",
                "Onco_PredR_AvgCount": str(round(df_onco[df_onco["RvN"] > 0.5]["Cosmic"].mean(), 2)),

                "TSG_PredA_Vars": str(len(df_tsg[df_tsg['AD_prediction'].str.contains('Activating')].index)),
                "TSG_PredA_%": f"{round(len(df_tsg[df_tsg['AD_prediction'].str.contains('Activating')].index) / len(df_tsg.index) * 100, 2)}%",
                "TSG_PredA_AvgCount": str(round(df_tsg[df_tsg['AD_prediction'].str.contains('Activating')]["Cosmic"].mean(), 2)),
                "TSG_PredD_Vars": str(len(df_tsg[df_tsg['AD_prediction'].str.contains('Deactivating')].index)),
                "TSG_PredD_%": f"{round(len(df_tsg[df_tsg['AD_prediction'].str.contains('Deactivating')].index) / len(df_tsg.index) * 100, 2)}%",
                "TSG_PredD_AvgCount": str(round(df_tsg[df_tsg['AD_prediction'].str.contains('Deactivating')]["Cosmic"].mean(), 2)),
                "TSG_PredR_Vars": str(len(df_tsg[df_tsg["RvN"] > 0.5].index)),
                "TSG_PredR_%": f"{round(len(df_tsg[df_tsg['RvN'] > 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                "TSG_PredR_AvgCount": str(round(df_tsg[df_tsg["RvN"] > 0.5]["Cosmic"].mean(), 2))
            }
            final_df.append(d)

    final_df = pd.DataFrame(final_df)
    final_df.to_csv(file_path["stats"], sep="\t", index=False)

    columns = {
        "VarSet": "Variant data set: all variants ('all') or genome-wide screen variants only ('gws')",
        "MinCount": "Minimum number of times a variant must be observed (sample count) in COSMIC to be included in the analysis",
        "Vars": "Number of unique variants",
        "TotVars": "Total number of variants (sum of sample counts)",
        "Genes": "Number of genes",
        "Oncog": "Number of oncogenes",
        "TSG": "Number of tumor suppressor genes",
        "Both": "Number of genes that are both oncogenes and tumor suppressor genes",
        "Unknown": "Number of genes without known role in cancer",

        "Oncog_Vars": "Number of unique variants in oncogenes",
        "Oncog_%": "Fraction from total number of unique variants",
        "Oncog_AvgCount": "Average sample count of variants in oncogenes",
        "TSG_Vars": "Number of unique variants in tumor suppressor genes",
        "TSG_%": "Fraction from total number of unique variants",
        "TSG_AvgCount": "Average sample count of variants in tumor suppressor genes",

        "KnownA_Vars": "Number of unique known activating variants",
        "KnownA_%": "Fraction from total number of unique variants",
        "KnownA_TotCount": "Total sample count of known activating variants",
        "KnownA_TotCount_%": "Fraction from total variant count",
        "KnownA_AvgCount": "Average sample count of known activating variants",
        "KnownA_Genes": "Number of genes with known activating variants",
        "KnownD_Vars": "Number of unique known deactivating variants",
        "KnownD_%": "Fraction from total number of unique variants",
        "KnownD_TotCount": "Total sample count of known deactivating variants",
        "KnownD_TotCount_%": "Fraction from total variant count",
        "KnownD_AvgCount": "Average sample count of known deactivating variants",
        "KnownD_Genes": "Number of genes with known deactivating variants",
        "KnownR_Vars": "Number of unique known resistance variants",
        "KnownR_%": "Fraction from total number of unique variants",
        "KnownR_TotCount": "Total sample count of known resistance variants",
        "KnownR_TotCount_%": "Fraction from total variant count",
        "KnownR_AvgCount": "Average sample count of known resistance variants",
        "KnownR_Genes": "Number of genes with known resistance variants",
        "KnownN_Vars": "Number of unique known neutral variants",
        "KnownN_%": "Fraction from total number of unique variants",
        "KnownN_TotCount": "Total sample count of known neutral variants",
        "KnownN_TotCount_%": "Fraction from total variant count",
        "KnownN_AvgCount": "Average sample count of known neutral variants",
        "KnownN_Genes": "Number of genes with known neutral variants",

        "PredA_Vars": "Number of unique variants predicted as activating",
        "PredA_%": "Fraction from total number of unique variants",
        "PredA_TotCount": "Total sample count of variants predicted as activating",
        "PredA_TotCount_%": "Fraction from total variant count",
        "PredA_AvgCount": "Average sample count of variants predicted as activating",
        "PredA_Genes": "Number of genes with variants predicted as activating",

        "PredD_Vars": "Number of unique variants predicted as deactivating",
        "PredD_%": "Fraction from total number of unique variants",
        "PredD_TotCount": "Total sample count of variants predicted as deactivating",
        "PredD_TotCount_%": "Fraction from total variant count",
        "PredD_AvgCount": "Average sample count of variants predicted as deactivating",
        "PredD_Genes": "Number of genes with variants predicted as deactivating",

        "PredU_Vars": "Number of unique variants predicted as 'uncertain'",
        "PredU_%": "Fraction from total number of unique variants",
        "PredU_TotCount": "Total sample count of variants predicted as 'uncertain'",
        "PredU_TotCount_%": "Fraction from total variant count",
        "PredU_AvgCount": "Average sample count of variants predicted as 'uncertain'",
        "PredU_Genes": "Number of genes with variants predicted as 'uncertain'",

        "PredR_Vars": "Number of unique variants predicted as resistance",
        "PredR_%": "Fraction from total number of unique variants",
        "PredR_TotCount": "Total sample count of variants predicted as resistance",
        "PredR_TotCount_%": "Fraction from total variant count",
        "PredR_AvgCount": "Average sample count of variants predicted as resistance",
        "PredR_Genes": "Number of genes with variants predicted as resistance",

        "Onco_PredA_Vars": "Number of unique variants predicted as activating in oncogenes",
        "Onco_PredA_%": "Fraction from total number of unique variants",
        "Onco_PredA_AvgCount": "Average sample count of variants predicted as activating in oncogenes",
        "Onco_PredD_Vars": "Number of unique variants predicted as deactivating in oncogenes",
        "Onco_PredD_%": "Fraction from total number of unique variants",
        "Onco_PredD_AvgCount": "Average sample count of variants predicted as deactivating in oncogenes",
        "Onco_PredR_Vars": "Number of unique variants predicted as resistance in oncogenes",
        "Onco_PredR_%": "Fraction from total number of unique variants",
        "Onco_PredR_AvgCount": "Average sample count of variants predicted as resistance in oncogenes",
        "TSG_PredA_Vars": "Number of unique variants predicted as activating in tumor suppressor genes",
        "TSG_PredA_%": "Fraction from total number of unique variants",
        "TSG_PredA_AvgCount": "Average sample count of variants predicted as activating in tumor suppressor genes",
        "TSG_PredD_Vars": "Number of unique variants predicted as deactivating in tumor suppressor genes",
        "TSG_PredD_%": "Fraction from total number of unique variants",
        "TSG_PredD_AvgCount": "Average sample count of variants predicted as deactivating in tumor suppressor genes",
        "TSG_PredR_Vars": "Number of unique variants predicted as resistance in tumor suppressor genes",
        "TSG_PredR_%": "Fraction from total number of unique variants",
        "TSG_PredR_AvgCount": "Average sample count of variants predicted as resistance in tumor suppressor genes"
    }
    with open(file_path["cols"], "w") as f:
        for k, v in columns.items():
            f.write(f"{k}\t{v}\n")


if __name__ == '__main__':
    file_path = get_data_paths()

    compute_stats_from_cosmic_ml(file_path)

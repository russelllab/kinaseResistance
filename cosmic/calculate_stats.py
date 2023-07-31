import sys
import pandas as pd
import numpy as np


def get_data_paths():
    """
    Get paths to all data files
    """
    file_paths = {
        "cosmic_ml_all": "ML_output_cosmic_all.tsv.gz",
        "cosmic_ml_gws": "ML_output_cosmic_gws.tsv.gz",
        "stats": "stats_from_cosmic_activark.txt"
    }
    return file_paths


def print_stats_from_cosmic_ml(file_path):
    """
    Calculate a bunch of stats from the cosmic ML output files. Feel free to add more to it.
    Output is a df printed in a tab-delimited file.
    """
    final_df = []
    for mode, file in zip(["all", "gws"], [file_path["cosmic_ml_all"], file_path["cosmic_ml_gws"]]):
        df = pd.read_csv(file, sep="\t")
        df['GeneRole'] = df['GeneRole'].replace(np.nan, "-")

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
            pred_act = df2[(df2["A"] > df2["D"]) & (df2["A"] > df2["N"])]
            pred_deact = df2[(df2["D"] > df2["A"]) & (df2["D"] > df2["N"])]
            pred_neut = df2[(df2["N"] > df2["A"]) & (df2["N"] > df2["D"])]
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

                "PredN_Vars": str(len(pred_neut.index)),
                "PredN_%": f"{round(len(pred_neut.index) / len(df2.index) * 100, 2)}%",
                "PredN_TotCount": str(round(pred_neut["Cosmic"].sum(), 2)),
                "PredN_TotCount_%": f"{round(pred_neut['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "PredN_AvgCount": str(round(pred_neut["Cosmic"].mean(), 2)),
                "PredN_Genes": str(pred_neut["GeneName"].nunique()),

                "PredR_Vars": str(len(pred_res.index)),
                "PredR_%": f"{round(len(pred_res.index) / len(df2.index) * 100, 2)}%",
                "PredR_TotCount": str(round(pred_res["Cosmic"].sum(), 2)),
                "PredR_TotCount_%": f"{round(pred_res['Cosmic'].sum() / df2['Cosmic'].sum() * 100, 2)}%",
                "PredR_AvgCount": str(round(pred_res["Cosmic"].mean(), 2)),
                "PredR_Genes": str(pred_res["GeneName"].nunique()),

                "Res_Act_Vars": str(len(pred_res[(pred_res["A"] > pred_res["D"]) & (pred_res["A"] > pred_res["N"])])),
                "Res_Act_%": f"{round(len(pred_res[(pred_res['A'] > pred_res['D']) & (pred_res['A'] > pred_res['N'])]) / len(df2.index) * 100, 2)}%",
                "Res_Deact_Vars": str(len(pred_res[(pred_res["D"] > pred_res["A"]) & (pred_res["D"] > pred_res["N"])])),
                "Res_Deact_%": f"{round(len(pred_res[(pred_res['D'] > pred_res['A']) & (pred_res['D'] > pred_res['N'])]) / len(df2.index) * 100, 2)}%",
                "Res_Neut_Vars": str(len(pred_res[(pred_res["N"] > pred_res["A"]) & (pred_res["N"] > pred_res["D"])])),
                "Res_Neut_%": f"{round(len(pred_res[(pred_res['N'] > pred_res['A']) & (pred_res['N'] > pred_res['D'])]) / len(df2.index) * 100, 2)}%",

                # "Onco_PredA_Vars": str(len(df_onco[df_onco["AIvLD"] > 0.5].index)),
                # "Onco_PredA_%": f"{round(len(df_onco[df_onco['AIvLD'] > 0.5].index) / len(df_onco.index) * 100, 2)}%",
                # "Onco_PredA_AvgCount": str(round(df_onco[df_onco["AIvLD"] > 0.5]["Cosmic"].mean(), 2)),
                # "Onco_PredD_Vars": str(len(df_onco[df_onco["AIvLD"] <= 0.5].index)),
                # "Onco_PredD_%": f"{round(len(df_onco[df_onco['AIvLD'] <= 0.5].index) / len(df_onco.index) * 100, 2)}%",
                # "Onco_PredD_AvgCount": str(round(df_onco[df_onco["AIvLD"] <= 0.5]["Cosmic"].mean(), 2)),
                # "Onco_PredR_Vars": str(len(df_onco[df_onco["RvN"] > 0.5].index)),
                # "Onco_PredR_%": f"{round(len(df_onco[df_onco['RvN'] > 0.5].index) / len(df_onco.index) * 100, 2)}%",
                # "Onco_PredR_AvgCount": str(round(df_onco[df_onco["RvN"] > 0.5]["Cosmic"].mean(), 2)),

                # "TSG_PredA_Vars": str(len(df_tsg[df_tsg["AIvLD"] > 0.5].index)),
                # "TSG_PredA_%": f"{round(len(df_tsg[df_tsg['AIvLD'] > 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                # "TSG_PredA_AvgCount": str(round(df_tsg[df_tsg["AIvLD"] > 0.5]["Cosmic"].mean(), 2)),
                # "TSG_PredD_Vars": str(len(df_tsg[df_tsg["AIvLD"] <= 0.5].index)),
                # "TSG_PredD_%": f"{round(len(df_tsg[df_tsg['AIvLD'] <= 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                # "TSG_PredD_AvgCount": str(round(df_tsg[df_tsg["AIvLD"] <= 0.5]["Cosmic"].mean(), 2)),
                # "TSG_PredR_Vars": str(len(df_tsg[df_tsg["RvN"] > 0.5].index)),
                # "TSG_PredR_%": f"{round(len(df_tsg[df_tsg['RvN'] > 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                # "TSG_PredR_AvgCount": str(round(df_tsg[df_tsg["RvN"] > 0.5]["Cosmic"].mean(), 2))

                # "Top5Genes(NumVars)": ";".join(
                #     df2.groupby("GeneName").count()["Cosmic"].sort_values(ascending=False).head(5).index.values),
                # "Top5Genes(VarCount)": ";".join(
                #     df2.groupby("GeneName").sum()["Cosmic"].sort_values(ascending=False).head(5).index.values)
            }
            final_df.append(d)

    final_df = pd.DataFrame(final_df)
    final_df.to_csv(file_path["stats"], sep="\t", index=False)

if __name__ == '__main__':
    file_path = get_data_paths()

    print_stats_from_cosmic_ml(file_path)

import pandas as pd
import numpy as np

def get_data_paths():
    """
    Get paths to all data files
    Modify this function to point to the correct paths on your system
    """
    cosmic_path = "/net/home.isilon/ds-russell/COSMIC/latest_version/"
    uniprot_path = "/net/home.isilon/ds-russell/uniprot/knowledgebase/complete/"
    file_path = {
        "cosmic_gws": f"{cosmic_path}/CosmicGenomeScreensMutantExport_mech_gn_mapped.txt",
        "cosmic_tar": f"{cosmic_path}/CosmicCompleteTargetedScreensMutantExport_mech_gn_mapped.txt",
        "cosmic_census": f"{cosmic_path}/cancer_gene_census.csv",
        "uni_fasta_file": f"{uniprot_path}/uniprot_sprot.fasta.gz",
        "ml_output_dir": "/net/home.isilon/ds-russell/kinaseResistance/ML/outputs_old/",
        "ml_output_file": "cosmic_activark.txt.gz",
        "mut_counts_all": "mutation_counts_all.tsv.gz",
        "mut_counts_gws": "mutation_counts_gws.tsv.gz",
        "cosmic_ml_all": "ML_output_cosmic_all.tsv.gz",
        "cosmic_ml_gws": "ML_output_cosmic_gws.tsv.gz",
        "stats": "stats_from_cosmic_activark.txt"
    }
    return file_path

def print_stats_from_cosmicML(file_path):
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
            d = {
                "VarSet": mode,
                "MinCount": f">{count_min}",
                "Vars": str(len(df2.index)),
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

                "KnownA_Vars": str(len(df2[df2["KnownADR"].str.contains("activating") | df2["KnownADR"].str.contains("increase")].index)),
                "KnownA_%": f"{round(len(df2[df2['KnownADR'].str.contains('activating') | df2['KnownADR'].str.contains('increase')].index) / len(df2.index) * 100, 2)}%",
                "KnownA_AvgCount": str(round(df2[df2["KnownADR"].str.contains("activating") | df2["KnownADR"].str.contains("increase")]["Cosmic"].mean(), 2)),
               
                "KnownR_Vars": str(len(df2[df2["KnownADR"].str.contains("resistance")].index)),
                "KnownR_%": f"{round(len(df2[df2['KnownADR'].str.contains('resistance')].index) / len(df2.index) * 100, 2)}%",
                "KnownR_AvgCount": str(round(df2[df2["KnownADR"].str.contains("resistance")]["Cosmic"].mean(), 2)),

                "KnownD_Vars": str(len(df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")].index)),
                "KnownD_%": f"{round(len(df2[df2['KnownADR'].str.contains('decrease') | df2['KnownADR'].str.contains('loss')].index) / len(df2.index) * 100, 2)}%",
                "KnownD_AvgCount": str(round(df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")]["Cosmic"].mean(), 2)),
                
                "KnownN_Vars": str(len(df2[df2["KnownADR"].str.contains("neutral")].index)),
                "KnownN_%": f"{round(len(df2[df2['KnownADR'].str.contains('neutral')].index) / len(df2.index) * 100, 2)}%",
                "KnownN_AvgCount": str(round(df2[df2["KnownADR"].str.contains("neutral")]["Cosmic"].mean(), 2)),

                "PredA_Vars": str(len(df2[df2["AIvLD"] > 0.5].index)),
                "PredA_%": f"{round(len(df2[df2['AIvLD'] > 0.5].index) / len(df2.index) * 100, 2)}%",
                "PredA_AvgCount": str(round(df2[df2["AIvLD"] > 0.5]["Cosmic"].mean(), 2)),
                "PredD_Vars": str(len(df2[df2["AIvLD"] <= 0.5].index)),
                "PredD_%": f"{round(len(df2[df2['AIvLD'] <= 0.5].index) / len(df2.index) * 100, 2)}%",
                "PredD_AvgCount": str(round(df2[df2["AIvLD"] <= 0.5]["Cosmic"].mean(), 2)),
                "PredR_Vars": str(len(df2[df2["RvN"] > 0.5].index)),
                "PredR_%": f"{round(len(df2[df2['RvN'] > 0.5].index) / len(df2.index) * 100, 2)}%",
                "PredR_AvgCount": str(round(df2[df2["RvN"] > 0.5]["Cosmic"].mean(), 2)),
  
                "Onco_PredA_Vars": str(len(df_onco[df_onco["AIvLD"] > 0.5].index)),
                "Onco_PredA_%": f"{round(len(df_onco[df_onco['AIvLD'] > 0.5].index) / len(df_onco.index) * 100, 2)}%",
                "Onco_PredA_AvgCount": str(round(df_onco[df_onco["AIvLD"] > 0.5]["Cosmic"].mean(), 2)),
                "Onco_PredD_Vars": str(len(df_onco[df_onco["AIvLD"] <= 0.5].index)),
                "Onco_PredD_%": f"{round(len(df_onco[df_onco['AIvLD'] <= 0.5].index) / len(df_onco.index) * 100, 2)}%",
                "Onco_PredD_AvgCount": str(round(df_onco[df_onco["AIvLD"] <= 0.5]["Cosmic"].mean(), 2)),
                "Onco_PredR_Vars": str(len(df_onco[df_onco["RvN"] > 0.5].index)),
                "Onco_PredR_%": f"{round(len(df_onco[df_onco['RvN'] > 0.5].index) / len(df_onco.index) * 100, 2)}%",
                "Onco_PredR_AvgCount": str(round(df_onco[df_onco["RvN"] > 0.5]["Cosmic"].mean(), 2)),

                "TSG_PredA_Vars": str(len(df_tsg[df_tsg["AIvLD"] > 0.5].index)),
                "TSG_PredA_%": f"{round(len(df_tsg[df_tsg['AIvLD'] > 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                "TSG_PredA_AvgCount": str(round(df_tsg[df_tsg["AIvLD"] > 0.5]["Cosmic"].mean(), 2)),
                "TSG_PredD_Vars": str(len(df_tsg[df_tsg["AIvLD"] <= 0.5].index)),
                "TSG_PredD_%": f"{round(len(df_tsg[df_tsg['AIvLD'] <= 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                "TSG_PredD_AvgCount": str(round(df_tsg[df_tsg["AIvLD"] <= 0.5]["Cosmic"].mean(), 2)),
                "TSG_PredR_Vars": str(len(df_tsg[df_tsg["RvN"] > 0.5].index)),
                "TSG_PredR_%": f"{round(len(df_tsg[df_tsg['RvN'] > 0.5].index) / len(df_tsg.index) * 100, 2)}%",
                "TSG_PredR_AvgCount": str(round(df_tsg[df_tsg["RvN"] > 0.5]["Cosmic"].mean(), 2))

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
    print_stats_from_cosmicML(file_path)


import os
import re
import sys
import glob
import gzip
import pandas as pd
from collections import defaultdict

def get_data_paths():
    cosmic_path = "/net/home.isilon/ds-russell/COSMIC/latest_version/"
    file_path = {
        "ml_output_dir": "/net/home.isilon/ds-russell/kinaseResistance/ML/outputs_old/",
        "cosmic_gws": f"{cosmic_path}/CosmicGenomeScreensMutantExport_mech_gn_mapped.txt",
        "cosmic_tar": f"{cosmic_path}/CosmicCompleteTargetedScreensMutantExport_mech_gn_mapped.txt",
        "cosmic_census": f"{cosmic_path}/cancer_gene_census.csv",
        "mut_counts_all": "mutation_counts_all.tsv.gz",
        "mut_counts_gws": "mutation_counts_gws.tsv.gz",
        "cosmic_ml_all": "ML_output_cosmic_all.tsv.gz",
        "cosmic_ml_gws": "ML_output_cosmic_gws.tsv.gz",
        "ml_output_file": "../ML/cosmic_ml_input_muts_activark_op.tsv.gz",
        "uni_fasta_file": "/net/home.isilon/ds-russell/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
    }
    return file_path

def get_list_of_kinases(path):
    kinases = []
    for file in glob.glob(path+"*"):
        accession = file.split("/")[-1].split(".")[0]
        kinases.append(accession)
    return kinases

def get_uni_seqs(file, kinases):
    uni_seqs = defaultdict(str)
    for line in gzip.open(file, "rt"):
        if line[0] == ">":
            accession = line.split("|")[1]
        else:
            if accession in kinases:
                uni_seqs[accession] += line.strip().replace(" ","")
    return uni_seqs

def parse_cosmic_mutations(path, kinases, uni_seqs, mech_samples=None):
    if mech_samples is None:
        mech_samples = defaultdict(set)
    s = set()
    for i, line in enumerate(open(path, "r")):
        if i == 0:
            continue
        t = line.split("\t")
        mech = t[0].split()[0]
        uni_ac = mech.split("/")[0]
        sample_id = t[5]
        n = len(t)

        status, change = None, None
        if "GenomeScreens" in path:
            change = t[21]
            status = t[27]
            if n == 37:
                change = t[22]
                status = t[28]
        elif "TargetedScreens" in path:
            change = t[21]
            status = t[28]
            if n == 38:
                change = t[22]
                status = t[29]
        s.add(status)
        if (uni_ac in kinases
                and status in ["Confirmed somatic variant", "Reported in another cancer sample as somatic"]
                and change == "Substitution - Missense"
                and mech != "ERROR"):
            wt, pos = re.search("(\w)(\d+)\w", mech.split("/")[1]).group(1, 2)
            if wt == uni_seqs[uni_ac][int(pos)-1]:
                mech_samples[mech].add(sample_id)
    print(s)
    return mech_samples

def parse_cosmic_census(census_file):
    df = pd.read_csv(census_file, sep=",", quotechar='"')

    return df

def create_mutation_counts_list(path, mech_samples):
    mech_counts = defaultdict(int)
    with gzip.open(path, "wt") as f:
        for mech in sorted(mech_samples):
            mech_counts[mech] = len(mech_samples[mech])
            f.write(f"{mech}\t{mech_counts[mech]}\n")
    return mech_counts

def get_cosmic_counts(paths, kinases):

    if os.path.isfile(paths["mut_counts_all"]) and os.path.isfile(paths["mut_counts_gws"]):
        with gzip.open(paths["mut_counts_all"], "rt") as f:
            mech_counts_all = {line.split("\t")[0]: int(line.split("\t")[1]) for line in f}
        with gzip.open(paths["mut_counts_gws"], "rt") as f:
            mech_counts_gws = {line.split("\t")[0]: int(line.split("\t")[1]) for line in f}
    else:
        uni_seqs = get_uni_seqs(paths["uni_fasta_file"], kinases)
        mech_samples = parse_cosmic_mutations(paths["cosmic_gws"], kinases, uni_seqs)
        sys.exit()
        mech_counts_gws = create_mutation_counts_list(paths["mut_counts_gws"], mech_samples)

        mech_samples = parse_cosmic_mutations(paths["cosmic_tar"], kinases, uni_seqs, mech_samples=mech_samples)
        mech_counts_all = create_mutation_counts_list(paths["mut_counts_all"], mech_samples)
    return mech_counts_all, mech_counts_gws

# def create_cosmic_ml_output(file_path, kin, cosmic_counts):
#     if os.path.isfile(file_path["cosmic_ml"]):
#         return
#     first = True
#     with gzip.open(file_path["cosmic_ml"], "wt") as f:
#         for kinase in kin:
#             for line in gzip.open(file_path["ml_output_dir"]+kinase+".txt.gz", "rt"):
#                 if line[0] == "#":
#                     header = line.strip()
#                 else:
#                     t = line.strip().split("\t")
#                     mech = t[0]
#                     if mech in cosmic_counts:
#                         if first:
#                             f.write(f"{header}\tCosmic\n")
#                             first = False
#                         f.write(f"{line.strip()}\t{cosmic_counts[mech]}\n")

def create_cosmic_ml_output(ml_output_file, cosmic_ml_output_file, cosmic_counts, census):
    to_print = {}
    s = set()
    for i, line in enumerate(gzip.open(ml_output_file, "rt")):
        t = [x.replace(" ","") for x in line.strip().replace("C-term Kinase-domain", "C-term-Kinase-domain").split()]
        if len(t) == 13:
            t = t[:8] + ["-"] + t[8:]

        if i == 0:
            header = "\t".join(t+["Cosmic", "GeneRole"])
            continue
        mech = t[0].replace(" ", "")
        gene = t[2].replace(" ", "")
        role = "NA"
        if gene in census["Gene Symbol"].values:
            role = census.loc[census["Gene Symbol"] == gene, "Role in Cancer"].values[0]
        if t[-1] != "NA":
            if mech in cosmic_counts:
                to_print[ "\t".join( t+[str(cosmic_counts[mech])]+[str(role)] ) ] = cosmic_counts[mech]

    with gzip.open(cosmic_ml_output_file, "wt") as f:
        f.write(f"{header}\n")
        for row in sorted(to_print, key=to_print.get, reverse=True):
            f.write(row+"\n")


if __name__ == '__main__':
    file_path = get_data_paths()
    kinase_list = get_list_of_kinases(file_path["ml_output_dir"])
    cosmic_counts_all, cosmic_counts_gws = get_cosmic_counts(file_path, kinase_list)
    census = parse_cosmic_census(file_path["cosmic_census"])

    rerun = False
    if os.path.isfile(file_path["ml_output_file"]) and rerun == True:
        create_cosmic_ml_output(file_path["ml_output_file"], file_path["cosmic_ml_all"], cosmic_counts_all, census)
        create_cosmic_ml_output(file_path["ml_output_file"], file_path["cosmic_ml_gws"], cosmic_counts_gws, census)

    #### FIX MAPPING ####
    final_df = []
    for mode, file in zip(["all", "gws"], [file_path["cosmic_ml_all"], file_path["cosmic_ml_gws"]]):
        df = pd.read_csv(file, sep="\t")
        #print(df.head())
        for count_min in [0, 1, 2, 5, 10, 20 , 30]:
            df2 = df[df['Cosmic'] > count_min]
            d = {
                "VarSet": mode,
                "MinCount": f">{count_min}",
                "Vars": str(len(df2.index)),
                "KnownA_Vars": str(len(df2[df2["KnownADR"].str.contains("activ")].index)),
                "KnownR_Vars": str(len(df2[df2["KnownADR"].str.contains("resistance")].index)),
                "KnownD_Vars": str(len(df2[df2["KnownADR"].str.contains("decrease") | df2["KnownADR"].str.contains("loss")].index)),
                "NumberGenes": str(df2["GeneName"].nunique()),
                "Top5Genes (NumVars)": ";".join(
                    df2.groupby("GeneName").count()["Cosmic"].sort_values(ascending=False).head(5).index.values),
                "Top5Genes (VarCount)": ";".join(
                    df2.groupby("GeneName").sum()["Cosmic"].sort_values(ascending=False).head(5).index.values)

            }
            final_df.append(d)

    final_df = pd.DataFrame(final_df)
    print(final_df)
    final_df.to_csv("stats.txt", sep="\t", index=False)


        # kin = set()
        # kin_positions = set()
        # for mech in sorted(cosmic_counts, key=cosmic_counts.get, reverse=True):
        #     kin.add(mech.split("/")[0])
        #     kin_positions.add(mech[:-1])
        #     if cosmic_counts[mech] < 2:
        #         break
        #     # print(f"{mech}\t{cosmic_counts[mech]}")
        # sys.exit()
        # print(f"# {len(cosmic_counts)} variants")
        # print(f"# in {len(kin)} kinases")
        # print(f"# in {len(kin_positions)} kinase positions")
        #
        # create_cosmic_ml_output(file_path, kin, cosmic_counts)
        #
        # # Open Ml output file and sort all rows by cosmic count
        # with gzip.open(file_path["cosmic_ml"], "rt") as f:
        #     header = f.readline().strip()
        #     lines = f.readlines()
        # lines.sort(key=lambda x: int(x.strip().split("\t")[-1]), reverse=True)
        # for line in lines:
        #     print(line.strip())


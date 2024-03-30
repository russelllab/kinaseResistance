import os
import csv
import gzip
import pandas as pd


def get_data_paths():
    """
    Get paths to all data files
    """
    file_paths = {
        # Census file as provided by COSMIC
        "cosmic_census": "data/Cosmic_v98_cancer_gene_census.csv.gz",
       
        # These 2 files need to be generated before using https://github.com/jcgonzs/cosmic_tools
        "cosmic_counts": "data/Cosmic_v98_counts.tsv.gz",
        "cosmic_miss": "data/Cosmic_v98_mismatches.tsv.gz",
       
        # File listing all kinase UniProt accessions
        "kinase_list": "../ML/all_kinases_acc.txt.gz",
        
        # Machine learning predictions
        "ml_output_file": "ML/cosmic_activark.txt.gz",
        
        # Final output files
        "cosmic_ml_all": "results/ML_output_cosmic_all.tsv.gz",
        "cosmic_ml_gws": "results/ML_output_cosmic_gws.tsv.gz"
    }
    return file_paths


def get_list_of_kinases(path):
    """
    Get list of UniProt accessions of all human kinases for which there are predictions.
    :param path: path to directory with predictions for each kinase
    """
    kinases = []
    for line in gzip.open(path, "rt"):
        kinases.append(line.rstrip())
    return kinases


def get_cosmic_counts(counts_file, miss_file, kinases):
    """
    Get sample counts of COSMIC mutations for all kinases. Required files are generating by parsing COSMIC mutant
    export files using https://github.com/jcgonzs/cosmic_tools
    :param counts_file: path to COSMIC counts file
    :param miss_file: path to COSMIC mismatches file
    :param kinases: list of UniProt accessions of kinases
    """
    mismatches = []
    for row in csv.reader(gzip.open(miss_file, "rt"), delimiter="\t"):
        mismatches.append(row[0])

    mech_counts_all = {}
    mech_counts_gws = {}
    for row in csv.reader(gzip.open(counts_file, "rt"), delimiter="\t"):
        acc = row[0].split("/")[0]
        ensp = row[1]
        if acc in kinases and ensp not in mismatches:
            mech_counts_all[row[0]] = int(row[4])
            mech_counts_gws[row[0]] = int(row[5])
    return mech_counts_all, mech_counts_gws


def parse_cosmic_census(census_file):
    """
    Read COSMIC cancer gene census file.
    :param census_file: path to COSMIC cancer gene census file
    """
    df = pd.read_csv(census_file, sep=",", quotechar='"')
    df["Role in Cancer"] = df["Role in Cancer"].fillna("inconclusive")
    return df


def merge_cosmic_and_ml(ml_output_file, cosmic_ml_output_file, cosmic_counts, census):
    """
    Merge ML predictions with COSMIC mutation counts and cancer gene status.
    :param ml_output_file: path to file where final output is printed
    :param cosmic_ml_output_file: path to file with ML predictions for COSMIC mutations
    :param cosmic_counts: dictionary of COSMIC mutation sample counts
    :param census: COSMIC cancer gene census dataframe
    """
    to_print = {}

    for i, line in enumerate(gzip.open(ml_output_file, "rt")):
        if line[0] == "#":
            continue
        line = line.strip().replace("C-term Kinase-domain", "C-term-Kinase-domain")
        t = [x.replace(" ", "") for x in line.split()]

        if "UserInput" in line:
            header = "\t".join(t + ["Cosmic", "GeneRole"])
            continue

        mech = t[0]
        gene = t[2]
        role = "NA"
        if gene in census["Gene Symbol"].values:
            role = census.loc[census["Gene Symbol"] == gene, "Role in Cancer"].values[0]
        if t[-1] != "NA":
            if mech in cosmic_counts:
                to_print["\t".join(t + [str(cosmic_counts[mech])] + [str(role)])] = cosmic_counts[mech]

    with gzip.open(cosmic_ml_output_file, "wt") as f:
        f.write(f"{header}\n")
        for row in sorted(to_print, key=to_print.get, reverse=True):
            f.write(row + "\n")


if __name__ == '__main__':
    file_path = get_data_paths()

    kinase_list = get_list_of_kinases(file_path["kinase_list"])

    cosmic_counts_all, cosmic_counts_gws = get_cosmic_counts(file_path["cosmic_counts"],
                                                             file_path["cosmic_miss"],
                                                             kinase_list)

    census_gene_list = parse_cosmic_census(file_path["cosmic_census"])

    rerun = True
    if os.path.isfile(file_path["ml_output_file"]) and rerun is True:
        merge_cosmic_and_ml(file_path["ml_output_file"], file_path["cosmic_ml_all"], cosmic_counts_all,
                            census_gene_list)
        merge_cosmic_and_ml(file_path["ml_output_file"], file_path["cosmic_ml_gws"], cosmic_counts_gws,
                            census_gene_list)

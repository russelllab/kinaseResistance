
import os
import re
import sys
import csv
import glob
import gzip
import pandas as pd
from collections import defaultdict
from Bio import SeqIO


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
        "cosmic_gws_new": f"{cosmic_path}/Cosmic_GenomeScreensMutant_v98_GRCh38.tsv.gz",
        "cosmic_tar_new": f"{cosmic_path}/Cosmic_CompleteTargetedScreensMutant_v98_GRCh38.tsv.gz",
        "cosmic_census": f"{cosmic_path}/cancer_gene_census.csv",
        "uni_fasta_file": f"{uniprot_path}/uniprot_sprot.fasta.gz",
        "ml_output_dir": "/net/home.isilon/ds-russell/kinaseResistance/ML/outputs_old/",
        "ml_output_file": "cosmic_activark.txt.gz",
        "mut_counts_all": "mutation_counts_all.tsv.gz",
        "mut_counts_gws": "mutation_counts_gws.tsv.gz",
        "mut_counts_all_new": "mutation_counts_all_v98.tsv.gz",
        "mut_counts_gws_new": "mutation_counts_gws_v98.tsv.gz",
        "cosmic_ml_all": "ML_output_cosmic_all.tsv.gz",
        "cosmic_ml_gws": "ML_output_cosmic_gws.tsv.gz"
    }
    return file_path


def get_list_of_kinases(path):
    """
    Get list of UniProt accessions of all human kinases for which there are
    predictions.
    """
    kinases = []
    for file in glob.glob(path+"*"):
        accession = file.split("/")[-1].split(".")[0]
        kinases.append(accession)
    return kinases


def get_cosmic_counts(paths, kinases, new=False):
    """
    Get sample counts of COSMIC mutations for all kinases. If counts files do 
    not exist, create them.
    :param paths: dictionary of paths to data files
    :param kinases: list of UniProt accessions of kinases
    :param new: whether to use new COSMIC files
    """
    if new:
        paths["cosmic_gws"] = paths["cosmic_gws_new"]
        paths["cosmic_tar"] = paths["cosmic_tar_new"]
        paths["mut_counts_all"] = paths["mut_counts_all_new"]
        paths["mut_counts_gws"] = paths["mut_counts_gws_new"]
    if os.path.isfile(paths["mut_counts_all"]) and os.path.isfile(paths["mut_counts_gws"]):
        with gzip.open(paths["mut_counts_all"], "rt") as f:
            mech_counts_all = {line.split("\t")[0]: int(line.split("\t")[1]) for line in f}
        with gzip.open(paths["mut_counts_gws"], "rt") as f:
            mech_counts_gws = {line.split("\t")[0]: int(line.split("\t")[1]) for line in f}
    else:
        uni_seqs = get_uni_seqs(paths["uni_fasta_file"], kinases)
        mech_samples = parse_cosmic_mutations(paths["cosmic_gws"], kinases, uni_seqs)
        mech_counts_gws = create_mutation_counts_list(paths["mut_counts_gws"], mech_samples)
        mech_samples = parse_cosmic_mutations(paths["cosmic_tar"], kinases, uni_seqs, mech_samples=mech_samples)
        mech_counts_all = create_mutation_counts_list(paths["mut_counts_all"], mech_samples)
    return mech_counts_all, mech_counts_gws


def get_uni_seqs(file, kinases=None):
    """
    Get UniProt sequences. If provided, restrict to provided list of accessions.
    :param file: path to UniProt FASTA file
    :param kinases: list of UniProt accessions of kinases
    """
    uni_seqs = {}
    for record in SeqIO.parse(gzip.open(file, "rt"), "fasta"):
        accession = record.id.split("|")[1]
        if kinases is None or accession in kinases:
            uni_seqs[accession] = str(record.seq)
    return uni_seqs


def parse_cosmic_mutations(path, kinases, uni_seqs, mech_samples=None):
    """
    Parse somatic missense variants from COSMIC mutant export file.
    If dictionary of mutation samples is provided, add to it.
    :param path: path to COSMIC mutant export file
    :param kinases: list of UniProt accessions of kinases
    :param uni_seqs: dictionary of UniProt sequences
    :param mech_samples: dictionary of COSMIC mutation counts
    """   
    if mech_samples is None:
        mech_samples = defaultdict(set)

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

        if (uni_ac in kinases
        and status in ["Confirmed somatic variant", "Reported in another cancer sample as somatic"]
        and change == "Substitution - Missense"
        and mech != "ERROR"):
            wt, pos = re.search("(\w)(\d+)\w", mech.split("/")[1]).group(1, 2)
            if wt == uni_seqs[uni_ac][int(pos)-1]:
                mech_samples[mech].add(sample_id)

    return mech_samples


def create_mutation_counts_list(path, mech_samples):
    """
    Write file with COSMIC mutation sample counts for each kinase and return 
    dictionary of counts.
    :param path: path to output file
    :param mech_samples: dictionary of COSMIC mutation and list of samples
    :return: dictionary of COSMIC mutation sample counts
    """
    mech_counts = defaultdict(int)
    with gzip.open(path, "wt") as f:
        for mech in sorted(mech_samples):
            mech_counts[mech] = len(mech_samples[mech])
            f.write(f"{mech}\t{mech_counts[mech]}\n")
    return mech_counts


def parse_cosmic_census(census_file):
    """
    Read COSMIC cancer gene census file.
    :param census_file: path to COSMIC cancer gene census file
    """
    df = pd.read_csv(census_file, sep=",", quotechar='"')

    return df


def merge_cosmic_MLoutput(ml_output_file, cosmic_ml_output_file, cosmic_counts, 
                          census):
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
        t = [x.replace(" ","") for x in line.split()]
        
        if "UserInput" in line:
            header = "\t".join(t+["Cosmic", "GeneRole"])
            continue
        
        mech = t[0]
        gene = t[2]
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
    cosmic_counts_all, cosmic_counts_gws = get_cosmic_counts(file_path, 
                                                            kinase_list, new=True)
    census = parse_cosmic_census(file_path["cosmic_census"])

    rerun = False
    if os.path.isfile(file_path["ml_output_file"]) and rerun == True:
        merge_cosmic_MLoutput(file_path["ml_output_file"], file_path["cosmic_ml_all"], 
                              cosmic_counts_all, census)
        merge_cosmic_MLoutput(file_path["ml_output_file"], file_path["cosmic_ml_gws"], 
                              cosmic_counts_gws, census)

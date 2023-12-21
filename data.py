from typing import Any
from Bio import Entrez
from Bio import Align
import numpy as np
import pandas as pd
from requests import get, codes as http_code
from graph import Graph
from mutation import Mutation
import networkx as nx
from pyvis.network import Network


def clean_data(df):
    arr = df.to_numpy()
    num_vars = len(arr)
    drop_rows = set()

    # col 4 describes pathogenicity
    for i in range(num_vars):
        curr = arr[i][4]
        if "Pathogen" not in arr[i][4]:
            drop_rows.add(i)
        elif arr[i][13] is np.nan:
            drop_rows.add(i)
        elif "Cystic fibrosis" not in arr[i][3]:
            drop_rows.add(i)
        elif "no assertion criteria provided" in arr[i][5]:
            drop_rows.add(i)
        if type(arr[i][13]) is str:
            if '|' in arr[i][13]:
                drop_rows.add(i)
                '''
                rs = ""
                j = 0
                s = arr[i][13][j]
                while s != "|":
                    rs += s
                    j += 1
                    s = arr[i][13][j]

                arr[i][13] = rs '''

        if type(arr[i][14]) is float:
            drop_rows.add(i)

    drop_rows = list(drop_rows)
    df = df.drop(labels=drop_rows, axis=0)
    arr = df.to_numpy()
    return arr


def query_by_rsid(rsid):
    # this function is adapted from the ncbi's dbsnp
    # purpose is to acquire json object of the SNP using the rsid
    access = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + str(rsid[2:]) + "/frequency"
    reply = get(access)
    if reply.status_code != http_code.ok:
        return

    content_type = reply.headers['content-type']
    if content_type != 'application/json':
        raise Exception("Unexpected content type: {}\n{}".format(
            content_type, access))

    return reply.json()


def get_sequence_by_rsid(rsid):
    Entrez.email = "sereenabhanji20@gmail.com"

    handle = Entrez.esearch(db="snp", term=rsid, retmode="xml")
    record = Entrez.read(handle)

    if record['Count'] == '0':
        print("No records found for rsID {rsid}")
        return None

    snp_id = str(record['IdList'][0])
    handle = Entrez.efetch(db="snp", id=snp_id, retmode="xml")
    snp_record = Entrez.read(handle, validate=False)
    print(handle.readline().strip())

    try:
        sequence = snp_record[0]['Assay']['Sequence']
        return sequence
    except KeyError:
        print("No sequence information found for rsID {rsid}")
        return None


def print_study_counts(study_id: str, study_counts: Any) -> None:

    print("\tAllele counts for study: {}".format(study_id))
    allele_counts = study_counts["allele_counts"]

    for pop_id, pop_counts in allele_counts.items():
        print("\t\tAllele counts for population {}".format(pop_id))
        for allele, count in pop_counts.items():
            print("\t\t\tAllele: {}. Count: {}".format(allele, count))


def get_snp_ranges(arr):
    ranges = []
    for snp in arr:
        loc = snp[10]
        if "-" in loc:
            start_pos = ""
            index = 0
            curr = loc[index]
            while curr != " ":
                start_pos += curr
                index += 1
                curr = loc[index]
            start_pos = int(start_pos)

            end_pos = int(loc[(index + 2):len(loc)])
        else:
            start_pos = int(loc)
            end_pos = int(loc)

        sequence = snp[14]
        index = 0
        while sequence[len(sequence) - 1 - index] != ":":
            index += 1
        alt_allele = sequence[(len(sequence) - index):]

        mut = Mutation(snp[13], start_pos, end_pos, alt_allele)
        ranges.append(mut)
    return ranges


def construct_graph(ranges):
    size = len(ranges)
    g = Graph(size)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    for i in range(size):
        for j in range(size):
            if i != j:
                m1 = ranges[i]
                m2 = ranges[j]
                weight = 0
                if m1.start <= m2.end & m2.start <= m1.end:
                    w = min(m1.end, m2.end) - max(m1.start, m2.start)
                    if w >= 0:
                        weight += w

                    s1 = m1.sequence
                    s2 = m2.sequence
                    if s1 != '' and s2 != '':
                        alignment = aligner.align(s1, s2)
                        score = alignment.score
                        weight += int(score)

                    if weight >= 1:
                        g.add_edge((i, j, weight))
    return g


def visualize_graph(g):
    FG = nx.Graph()
    for edge in g:
        FG.add_edge(edge[0], edge[1], weight=edge[2], label=str(edge[2]))

    nt = Network(notebook=True)
    nt.from_nx(FG)
    nt.show("chart3.html")


def mst_phylogeny(csv):

    df = pd.read_csv(csv, sep="\t")
    arr = clean_data(df)
    ranges = get_snp_ranges(arr)
    g = construct_graph(ranges)
    mst = g.kruskal()
    visualize_graph(mst)


if __name__ == '__main__':
    mst_phylogeny("/Users/sereenabhanji/Downloads/clinvar_result.txt")

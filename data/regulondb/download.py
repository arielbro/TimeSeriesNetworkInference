"""Download the RegulonDB E. coli K-12 regulatory network and the GSE20305 stress time series into
data/regulondb/raw, and write the network as an unsigned edge list with its signs in a separate file.
See README.md for what each file holds and how it was processed.

Sources:
    RegulonDB (regulondb.ccg.unam.mx), through its GraphQL service's getDataOfFile query:
        NetworkRegulatorGene           regulator -> regulated gene, with effect and confidence level
        GeneProductAllIdentifiersSet   gene names, synonyms and b-numbers
    GEO GSE20305 (Jozefczuk et al. 2010, Mol Syst Biol 6:364), from the NCBI FTP site:
        GSE20305_series_matrix.txt.gz  122 arrays, normalized intensities per probeset
        GPL10041 platform table        probeset -> ORF name and b-number, saved as a small tsv

Written to data/regulondb:
    network.tsv   regulator  target  regulator_class  confidence   (one row per gene pair, no signs)
    signs.tsv     regulator  target  sign   with sign 1, -1, dual or unknown
"""

import argparse
import collections
import gzip
import io
import json
import os
import urllib.request

REGULONDB_GRAPHQL = "https://regulondb.ccg.unam.mx/graphql"
REGULONDB_FILES = ["NetworkRegulatorGene", "GeneProductAllIdentifiersSet"]
GEO_MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20305/matrix/GSE20305_series_matrix.txt.gz"
GEO_PLATFORM = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10041/soft/GPL10041_family.soft.gz"
OUT_DIR = os.path.join("data", "regulondb")

# RegulonDB regulator id prefixes: TFC transcription factor, PDC other protein (mostly post-transcriptional,
# e.g. CsrA, Hfq, ribosomal proteins), RCC small RNA, CNC compound (only ppGpp, which has no gene and is dropped)
REGULATOR_CLASSES = {"TFC": "tf", "PDC": "protein", "RCC": "srna"}
SIGNS = {"+": "1", "-": "-1", "-+": "dual", "+-": "dual", "": "unknown", "?": "unknown"}
CONFIDENCE_ORDER = "?WSC"  # unknown < weak < strong < confirmed


def fetch(url, data=None):
    request = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"} if data else {})
    with urllib.request.urlopen(request, timeout=600) as response:
        return response.read()


def fetch_regulondb_file(name):
    """(file content, {version, creationDate}) for a RegulonDB downloadable file."""
    query = '{ getDataOfFile(fileName: "%s") { content rdbVersion creationDate } }' % name
    result = json.loads(fetch(REGULONDB_GRAPHQL, json.dumps({"query": query}).encode()))["data"]["getDataOfFile"]
    return result["content"], {"rdbVersion": result["rdbVersion"], "creationDate": result["creationDate"]}


def platform_table(soft_gz):
    """The GPL10041 platform rows as tsv text with the ID, ORF and synonym columns."""
    lines, inside = [], False
    for line in gzip.open(io.BytesIO(soft_gz), "rt", encoding="utf8", errors="replace"):
        if line.startswith("!platform_table_begin"):
            inside = True
        elif line.startswith("!platform_table_end"):
            break
        elif inside:
            fields = line.rstrip("\r\n").split("\t")
            lines.append("\t".join(fields[i] for i in (0, 2, 3)))
    return "\n".join(lines) + "\n"


def data_rows(content):
    """Tab-separated rows of a RegulonDB file, comment and column-header lines skipped."""
    return [[field.strip() for field in line.split("\t")] for line in content.splitlines()
            if line.strip() and not line.startswith("#") and not line.startswith("1)")]


def merge_signs(first, second):
    """The sign of a gene pair reported twice: an unknown sign defers to a known one, and two different known
    signs make it dual."""
    if first == second or second == "unknown":
        return first
    return second if first == "unknown" else "dual"


def build_network(content):
    """{(regulator gene, target gene): (class, confidence, sign)} at the gene level. A regulator coded by
    several genes (a complex such as IHF = ihfA;ihfB) gets one edge per subunit gene. When one gene pair
    arises from several RegulonDB rows (two complexes sharing a subunit), the highest confidence is kept
    and differing signs become dual."""
    edges = {}
    for regulator_id, _, regulator_genes, _, target, function, confidence in data_rows(content):
        regulator_class = REGULATOR_CLASSES.get(regulator_id[8:11])
        if regulator_class is None:
            continue
        sign = SIGNS[function]
        for regulator in [gene.strip() for gene in regulator_genes.replace(",", ";").split(";") if gene.strip()]:
            key = (regulator, target)
            if key in edges:
                old_class, old_confidence, old_sign = edges[key]
                sign = merge_signs(old_sign, sign)
                confidence = max(old_confidence, confidence, key=CONFIDENCE_ORDER.index)
            edges[key] = (regulator_class, confidence, sign)
    return edges


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=OUT_DIR)
    options = parser.parse_args()
    raw_dir = os.path.join(options.out, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    versions = {}
    for name in REGULONDB_FILES:
        content, versions[name] = fetch_regulondb_file(name)
        with open(os.path.join(raw_dir, name + ".tsv"), "w", encoding="utf8", newline="\n") as out_file:
            out_file.write(content)
    with open(os.path.join(raw_dir, "regulondb_versions.json"), "w") as out_file:
        json.dump(versions, out_file, indent=1)
    with open(os.path.join(raw_dir, "GSE20305_series_matrix.txt.gz"), "wb") as out_file:
        out_file.write(fetch(GEO_MATRIX))
    with open(os.path.join(raw_dir, "GPL10041_platform.tsv"), "w", encoding="utf8", newline="\n") as out_file:
        out_file.write(platform_table(fetch(GEO_PLATFORM)))

    with open(os.path.join(raw_dir, "NetworkRegulatorGene.tsv"), encoding="utf8") as in_file:
        edges = build_network(in_file.read())
    with open(os.path.join(options.out, "network.tsv"), "w", newline="\n") as network_file, \
            open(os.path.join(options.out, "signs.tsv"), "w", newline="\n") as signs_file:
        network_file.write("regulator\ttarget\tregulator_class\tconfidence\n")
        signs_file.write("regulator\ttarget\tsign\n")
        for (regulator, target), (regulator_class, confidence, sign) in sorted(edges.items()):
            network_file.write("\t".join([regulator, target, regulator_class, confidence]) + "\n")
            signs_file.write("\t".join([regulator, target, sign]) + "\n")

    genes = {gene for edge in edges for gene in edge}
    print("RegulonDB {}: {} edges over {} genes, {} regulators".format(
        versions["NetworkRegulatorGene"]["rdbVersion"], len(edges), len(genes), len({r for r, _ in edges})))
    print("regulator classes", dict(collections.Counter(c for c, _, _ in edges.values())))
    print("confidence", dict(collections.Counter(c for _, c, _ in edges.values())))
    print("signs", dict(collections.Counter(s for _, _, s in edges.values())))


if __name__ == "__main__":
    main()

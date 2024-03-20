import requests
import scipy.stats as stats
import pandas as pd
import numpy as np

def reformat_gs_to_https(gs_path, project_id):
    if gs_path.startswith('gs'):
        gs_path = gs_path[5:]
    elif not gs_path.startswith('fc'):
        raise ValueError(f'Provided path must begin with gs or fc: {gs_path}')
        
    auth_path = 'https://storage.cloud.google.com/' + gs_path + f'?userProject={project_id}'
    
    return auth_path


def read_gtf(gtf_url):
    
    gtf = pd.read_csv(
        gtf_url,
        comment="#",
        sep="\t",
        header=None,
    )

    gtf.columns = [
        "chromosome",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    gtf = gtf[gtf["feature"] == "gene"].drop_duplicates()
    gtf["attributes"] = gtf.attributes.apply(
        lambda x: {
            i.split(" ")[0]: i.split(" ")[1].replace('"', "")
            for i in x.strip(";").split("; ")
        }
    )
    gtf["gene_name"] = gtf.attributes.apply(lambda x: x["gene_name"])
    gtf["gene_id"] = gtf.attributes.apply(lambda x: x["gene_id"].split(".")[0])
    gtf["gene_type"] = gtf.attributes.apply(lambda x: x["gene_type"])
    
    return gtf


def clean_and_subset_gtf(gtf, gene_ref):
    '''Given a gene reference and a gtf, clean up the gtf'''
    
    cleaned_gtf = gtf.query('gene_type == "protein_coding"').copy()
    cleaned_gtf = cleaned_gtf[cleaned_gtf["gene_name"].isin(gene_ref)]

    # remove PAR genes
    par = np.array(["PAR" in row["attributes"]["gene_id"] for _, row in cleaned_gtf.iterrows()])
    cleaned_gtf = cleaned_gtf[~par]

    # handle duplicate genes. Keep level 1 entries and otherwise pick the first
    nonduplicated = cleaned_gtf[~cleaned_gtf["gene_name"].duplicated(keep=False)]
    duplicated = cleaned_gtf[cleaned_gtf["gene_name"].duplicated(keep=False)]

    level1 = np.array(
        [row["attributes"]["level"] == "1" for _, row in duplicated.iterrows()]
    )
    level1_genes = duplicated[level1]

    duplicated = duplicated[~(duplicated["gene_name"].isin(level1_genes["gene_name"]))]

    duplicated = duplicated.drop_duplicates(subset="gene_name", keep="last")

    cleaned_gtf = pd.concat([nonduplicated, level1_genes, duplicated])
    
    # get the index of these genes and indicate that they are kept
    kept_rows = sorted(cleaned_gtf.index)
    
    gtf['kept_for_subset'] = False
    gtf.loc[kept_rows, 'kept_for_subset'] = True
    
    return gtf


def determine_parabola(x, y):
    
    a, b, c = np.polyfit(x, y, 2)
    
    plot_x = np.linspace(x[0], x[2], 100)
    plot_y = a*plot_x**2 + b*plot_x + c
    return (plot_x, plot_y)


def fetch_track_overlap(track, coord):
    '''Given a UCSC track name, return a binary for if an element in that track overlaps that region.
       
       Params:
       -------
       track: str, name of UCSC track
       coord: tuple, (chrom, start, end)
       
       Returns:
       --------
       overlap: bool, True if track has an element that overlaps that coordinate
       '''
    
    if coord[0][:3] != "chr":
        raise ValueError('chrom must be prefixed with "chr"')
    
    chrom, start, end = coord
    
    url = f'https://api.genome.ucsc.edu/getData/track?genome=hg38;track={track};chrom={chrom};start={start};end={end}'
    resp = requests.get(url).json()
    
    if len(resp[track]):
        return 1
    else:
        return 0
        
    
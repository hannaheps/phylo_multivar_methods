from qiime2 import Artifact,Metadata
from qiime2.plugins.diversity.pipelines import beta
from qiime2.plugins.diversity.methods import pcoa
from os.path import exists,splitext
import pandas as pd
from skbio import OrdinationResults

def extract_pc_coords(pcoa):
    """Return a DataFrame of PC axes"""
    pass


if __name__ == "__main__":

    metric = "Aitchison"
    n_pcoa_dimensions = 3

    pcoa_filepath = "../output/raw_sim_microbiome_1000microbes_rarefied_1000_pcoa.qza"

    if not exists(pcoa_filepath):  
        print(f"PCoA file not found at {pcoa_filepath}. Generating...")  
        feature_table = Artifact.load("../output/raw_sim_microbiome_1000microbes_rarefied_1000.qza")
        beta_results = beta(table=feature_table, metric=metric.lower())
        beta_dm = beta_results.distance_matrix
        pcoa_results = pcoa(beta_dm,n_pcoa_dimensions)
        print(pcoa_results.pcoa)
        pcoa_axes = pcoa_results.pcoa
        pcoa_axes.save(pcoa_filepath)

    print(f"Loading PCoA results from {pcoa_filepath}")
    pcoa_axes = Artifact.load(pcoa_filepath)
    print(pcoa_axes)
    pcoa_ordination = pcoa_axes.view(OrdinationResults)
    pcoa_df = pcoa_ordination.samples.iloc[:,:]
    print(pcoa_df)
    print(type(pcoa_df))
    

    pcoa_df.to_csv(splitext(pcoa_filepath)[0]+".csv")

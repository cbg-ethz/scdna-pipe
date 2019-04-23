from secondary_analysis import SecondaryAnalysis
import argparse


# parser = argparse.ArgumentParser()

sa = SecondaryAnalysis(sample_name='MATIWAQ-T_scD_Ar1v1.3', output_path='MATIWAQ', h5_path='..//test_dir/cnv_data.h5', genes_path='../../required_files/genes_of_interest/melanoma_gene_coordinates.tsv')
sa.remove_tenx_genomics_artifacts(bins='/Users/mtuncel/git_repos/dna-pipeline/required_files/10x_artifacts/GRCh37_10x_artifacts_mask_whole_genome.tsv')
sa.apply_phenograph()
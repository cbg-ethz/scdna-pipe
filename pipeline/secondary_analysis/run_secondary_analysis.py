from secondary_analysis import SecondaryAnalysis
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample_name", required=True)
parser.add_argument("-o","--output_path", required=True)
parser.add_argument("-h5","--h5_path", required=True, help="path to the input h5")
parser.add_argument("-g","--genes_path", required=True, help="path to the genes of interest and their coordinates")
parser.add_argument("-b","--bins_to_remove", required=True, help="10x genomics technology artifact bins")


args = parser.parse_args()

sa = SecondaryAnalysis(sample_name=args.sample_name, output_path=args.output_path, h5_path=args.h5_path, genes_path=args.genes_path)
sa.remove_tenx_genomics_artifacts(bins=args.bins_to_remove)
sa.apply_phenograph()


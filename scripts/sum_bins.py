import argparse
import numpy as np
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mat", required=True, help="The matrix containing original bins")
parser.add_argument("-r", "--reg_sizes", required=True, help="Region sizes file")
parser.add_argument("-o","--output_path",required=False, default="./", help="path to the output")
parser.add_argument("-s", "--sample_name",required=False, default="", help="name of the sample")

args = parser.parse_args()

print("Reading the input...")
counts_per_bin = np.loadtxt(args.mat,delimiter='\t')
n_rows = counts_per_bin.shape[0]

region_sizes = np.loadtxt(args.reg_sizes, delimiter='\t', dtype=np.int32)
n_regions = len(region_sizes)

sum_region_sizes = np.sum(region_sizes,dtype=np.int32)

print("Summing bins...")
reg_avg_per_bin = np.zeros((n_rows, sum_region_sizes))

for i in tqdm(range(0, n_rows)):
    bin_count = 0
    for reg_size in region_sizes:
        region_sum = np.sum(counts_per_bin[i][bin_count: bin_count+reg_size])
        for j in range(bin_count, bin_count+reg_size):
            reg_avg_per_bin[i][j] = region_sum / reg_size
        bin_count += reg_size


print("writing output...")

np.savetxt(args.output_path + '/' + args.sample_name +"_segmented_counts.tsv", reg_avg_per_bin, delimiter='\t')

print("Summing up bins is done!")
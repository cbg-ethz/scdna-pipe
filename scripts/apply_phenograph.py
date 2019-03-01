import h5py
import argparse
import numpy as np
import pandas as pd
import phenograph
from collections import Counter
from sklearn.preprocessing import normalize

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="filtered counts input")
parser.add_argument("-o","--output_path",required=False, default="./", help="path to the output")
parser.add_argument("-s", "--sample_name",required=False, default="", help="name of the sample")

args = parser.parse_args()

# clustering/classification params
n_jobs = 16
# points in the knn neighbourhood are weighted by the distance
weight='distance'

filtered_counts = np.loadtxt(args.input)
normalized_filtered_counts = normalize(filtered_counts,axis=1, norm='l1')


n_cells = normalized_filtered_counts.shape[0]

print("n_cells: " + str(n_cells))
n_neighbours = int(n_cells/10)
print("n_neighbours: " + str(n_neighbours))
communities, graph, Q = phenograph.cluster(data=normalized_filtered_counts,k=n_neighbours,n_jobs=n_jobs, jaccard=True)

# computing the distance matrix from graph
arr = graph.toarray()
arr_full = arr+arr.T
np.fill_diagonal(arr_full, 1)
dist = (arr_full- arr_full.max())*(-1)
np.fill_diagonal(dist, 0)

print("shape of the distance matrix:")
print(dist.shape)

# write dist to file
# later use dist for all of the plots
dist_fname = args.output_path+'/'+args.sample_name+"_phenograph_distance.csv"
np.savetxt(fname=dist_fname, X=dist, delimiter=',')

print(communities) # one of the outputs

communities_df = pd.DataFrame(communities,columns=['cluster'])
communities_df['cell_barcode'] = communities_df.index
communities_df = communities_df[['cell_barcode','cluster']] # order the columns
communities_df.head()


communities_df.to_csv(args.output_path + '/' + args.sample_name + "_clusters_phenograph_assignment.tsv",sep='\t',index=False)

# write the modularity score, for stability
f = open( args.output_path + '/' + args.sample_name + "_clustering_score.txt", 'w' )
f.write(str(Q))
f.close()


cells_by_cluster = []

community_ids = sorted(list(Counter(communities)))

for cluster in community_ids:
    cells_by_cluster.append(filtered_counts[communities==cluster])

avg_clusters = [m.mean(0) for m in cells_by_cluster] # 2nd output

avg_clusters_df = pd.DataFrame(avg_clusters)

avg_clusters_df['cluster_ids'] = community_ids # add the community_ids

avg_clusters_df.to_csv(args.output_path + '/' + args.sample_name + "_clusters_phenograph_profiles.tsv",sep='\t',index=False, header=True)

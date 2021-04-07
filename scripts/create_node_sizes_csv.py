import argparse
import json
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "-t", "--tree_json_path", required=True,
)

args = parser.parse_args()

cluster_tree_json = args.tree_json_path
output_file = cluster_tree_json.split("__")[0] + "__cluster_tree_node_sizes.csv"

with open(cluster_tree_json) as json_file:
    cluster_tree = json.load(json_file)

node_labels = []
node_sizes = []
for node in cluster_tree:
    if cluster_tree[node]["size"] > 0:
        node_labels.append(cluster_tree[node]["label"])
        node_sizes.append(cluster_tree[node]["size"])
order = np.argsort(np.array(node_labels))
node_sizes = np.array(node_sizes)[order]

print("saving the node sizes...")
np.savetxt(output_file, node_sizes, delimiter=",", fmt="%d")

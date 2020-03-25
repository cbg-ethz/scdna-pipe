import argparse
import numpy as np
from scicone import *

parser = argparse.ArgumentParser()
parser.add_argument(
    "-t", "--input_tree", required=True, help="File containing learned tree"
)
parser.add_argument(
    "-d",
    "--input_segmented_data",
    required=True,
    help="File containing the segmented cell bin counts used to learn tree",
)
parser.add_argument(
    "-r",
    "--input_segmented_region_sizes",
    required=True,
    help="File containing the size of each region",
)
parser.add_argument(
    "-n",
    "--region_neutral_states_file",
    required=True,
    help="File containing the neutral state of each region",
)
parser.add_argument(
    "-b", "--scicone_binary", required=True, help="Path to SCICoNE binary"
)
parser.add_argument("-o", "--output_file", required=True, help="Output file")

args = parser.parse_args()

input_tree = args.input_tree
input_segmented_data = args.input_segmented_data
input_segmented_region_sizes = args.input_segmented_region_sizes
region_neutral_states_file = args.region_neutral_states_file
bin = args.scicone_binary
output_file = args.output_file

segmented_data = np.loadtxt(input_segmented_data, delimiter=",")
segmented_region_sizes = np.loadtxt(input_segmented_region_sizes, delimiter=",")
region_neutral_states = np.loadtxt(region_neutral_states_file, delimiter=",")

learned_tree = Tree(bin, "")
learned_tree.read_tree_str(input_tree)
fusion_tree = create_fusion_tree(learned_tree, region_neutral_states)
fusion_tree.learn_tree(
    segmented_data,
    segmented_region_sizes,
    n_iters=0,
    move_probs=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00],
    seed=42,
    postfix="fusion_tree",
    initial_tree=fusion_tree,
    region_neutral_states=region_neutral_states,
    copy_number_limit=10,
)

unique_assignments, unique_assignments_idx, tree_cluster_sizes = np.unique(
    fusion_tree.outsputs["cell_node_ids"], return_index=True, return_counts=True
)

# Map node ID to subclone ID
inferred_cnvs = learned_tree["inferred_cnvs"]

unique_cnvs, tree_node_sizes = np.unique(
    inferred_cnvs, axis=0, return_counts=True
)  # clone-wise profiles

# Sort clones by distance to diploid profile
dist_to_diploid = []
diploid_profile = np.ones([unique_cnvs.shape[1]]) * 2
for c_id in range(unique_cnvs.shape[0]):
    dist_to_diploid.append(np.linalg.norm(unique_cnvs[c_id] - diploid_profile))
order = np.argsort(dist_to_diploid)
unique_cnvs = unique_cnvs[order]
tree_node_sizes = tree_node_sizes[order]

labels = np.empty(inferred_cnvs.shape[0])
for c_id in range(unique_cnvs.shape[0]):
    cells = np.where(np.all(inferred_cnvs == unique_cnvs[c_id], axis=1))[0]
    labels[cells] = c_id

# Map node IDs to labels
_, idx = np.unique(labels)
node_label_map = dict(zip(original_cell_node_ids[idx], labels[idx]))

# Now add entries for the fusion nodes to the map by prefixing with 'F'
for node in list(node_label_map):
    node_label_map[str(int(node_label_map[node]) + 1000)] = "F" + node

n_cells = inferred_cnvs.shape[0]
new_cell_node_ids = np.empty(n_cells)
# Replace cell_node_ids from fusion tree with new labels
for node_id in np.unique(fusion_cell_node_ids):
    cells = np.where(np.all(fusion_cell_node_ids == node_id, axis=1))[0]
    new_cell_node_ids[cells] == node_label_map[node_id]

# Count number of cells attached to each node
uniques, counts = np.unique(new_cell_node_ids, return_counts=True)

# Print to file
with open(output_file, "w") as text_file:
    print(
        "# Number of cells assigned to each of the possible fusion clusters",
        file=text_file,
    )
    print(
        "Fusion clusters are labeled 'F\{subclone_id\}',"
        " where \{subclone_id\} indicates the ID of the subclones whose CNV profiles"
        " are plotted in the output PNG file \{sample_prefix\}__cluster_profile_overlapping.",
        file=text_file,
    )
    print(
        f"Only fusion clusters for which at least one cell attached to are shown. All the others"
        " are apparently not supported by the data.",
        file=text_file,
    )
    fcounts = 0
    for i, unique in enumerate(uniques):
        if unique.contains("F"):
            fcounts += counts[i]
            print(f"{unique}: {counts[i]}", file=text_file)
    print("", file=text_file)  # end with empty newline
    print(
        f"In total, there are {fcounts} possible doublets/fusion cells, out of {n_cells} ({fcounts/n_cells})."
    )

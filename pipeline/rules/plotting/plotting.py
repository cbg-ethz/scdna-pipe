from ..process_cnvs import *

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from tqdm import tqdm
import os

def tree_to_graphviz(tree_path, node_sizes=None, gene_labels=False, bin_gene_region_df=None, genes_to_highlight=None, highlight_color='red', max_genes_per_line=10, output_path=None):
    """
        reads the file containing trees converts it to graphviz format
        :param tree_path: path to the tree file.
        :param node_sizes: dictionary containing the size of each node
        :param gene_labels: whether to label nodes with genes
        :param bin_gene_region_df: Pandas DataFrame with gene-region correspondence
        :param genes_to_highlight: List containing genes to highlight
        :return: string object containing the graphviz formatted tree
    """
    with open(tree_path) as f:
        list_tree_file = list(f)

    graphviz_header = ["digraph { \n", "node [style=filled,color=\"#D4C0D6\"]"
                "edge [arrowhead=none, color=\"#602A86\"]"]

    graphviz_labels = []
    graphviz_links = []

    graphviz_labels.append("0[label=\"Neutral\"]") # root

    for line in list_tree_file:
        if line.startswith("node 0:"):
            continue
        elif line.startswith("node"):
            comma_splits = line.split(",")

            comma_first = re.split(" |:",comma_splits[0])
            node_id = comma_first[1]
            p_id = comma_first[4]
            comma_rest = comma_splits[1:]
            comma_rest[0] = comma_rest[0].lstrip('[')
            comma_rest[-1] = comma_rest[-1].rstrip(']\n')
            merged_labels = []
            [k_begin, previous_v] = (int(x) for x in comma_rest[0].split(":"))
            k_end = k_begin
            for term in comma_rest[1:]: # events vector
                [k,v] = (int(x) for x in term.split(":"))
                if k==k_end+1 and v==previous_v:
                    k_end = k # update the end
                else:
                    if k_begin == k_end:
                        merged_labels.append(f"{previous_v:+}R{k_begin}")
                    else:
                        merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")
                    k_begin = k_end = k
                previous_v = v
            # print the last one
            if k_begin == k_end:
                merged_labels.append(f"{previous_v:+}R{k_begin}")
            else:
                merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")

            str_merged_labels = " ".join(f"{x}<br/>" if i%10 == 0 and i>0 else str(x) for i, x in enumerate(merged_labels))
            if gene_labels and bin_gene_region_df is not None:
                node_str = " ".join(merged_labels) # "+1R75 +1R218:219 +1R221:223"
                str_merged_labels = convert_node_regions_to_genes(node_str, bin_gene_region_df,
                                    priority_only=True, genes_to_highlight=genes_to_highlight,
                                    highlight_color=highlight_color, max_genes_per_line=max_genes_per_line)
            # Add node size
            newline = "<br/>"
            if node_sizes is not None:
                node_size = node_sizes[node_id]
                str_merged_labels = str_merged_labels + " " + newline + " " + newline + str(int(node_size)) + " cell"
                if int(node_size) > 1:
                    str_merged_labels = str_merged_labels + "s"
                str_merged_labels = str_merged_labels + " "

            graphviz_labels.append(f"{node_id}[label=<{str_merged_labels}>]") # use < > to allow HTML
            graphviz_links.append(f"{p_id} -> {node_id}")

    txt = graphviz_header+graphviz_labels+graphviz_links+["}"]

    if output_path is not None:
        with open(output_path, "w") as file:
            for line in txt:
                file.write(f"{line}\n")

    return txt

def plot_tree_graphviz(tree_graphviz_path, output_path):
    try:
        format = output_path.split('.')[-1]
    except:
        format = 'png'

    try:
        cmd_output = subprocess.run(["dot", f"-T{format}:cairo", f"{tree_graphviz_path}", "-o", f"{output_path}"])
    except subprocess.SubprocessError as e:
        print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
    else:
        print(f"subprocess out: {cmd_output}")
        print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")

def plot_heatmap(gene_cn_df, is_imputed=None, output_path=None):
    if np.all(~np.isnan(gene_cn_df)):
        annot = np.array(gene_cn_df.astype(int).astype(str))
    else:
        annot = np.array(gene_cn_df.astype(str))
    annot[np.where(gene_cn_df==4)] = ['4+'] * len(np.where(gene_cn_df==4)[0])
    annot = annot.astype(str)

    figure_width = gene_cn_df.shape[0] / 2 + 1.5
    plt.figure(figsize=(8, figure_width))
    cmap = sns.color_palette("RdBu_r", 5)
    heatmap = sns.heatmap(
        gene_cn_df,
        annot=annot,
        cmap=cmap,
        vmin=0,
        vmax=4,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"ticks": [0, 1, 2, 3, 4]},
        fmt = ''
    )
    colorbar = heatmap.collections[0].colorbar
    colorbar.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
    colorbar.set_ticklabels(['0', '1', '2', '3', '4+'])
    heatmap.set_title("Copy number values of genes per subclone")
    heatmap.set_facecolor("#656565")
    plt.xlabel('Subclone')
    b, t = plt.ylim() # discover the values for bottom and top
    b += 0.5 # Add 0.5 to the bottom
    t -= 0.5 # Subtract 0.5 from the top
    plt.ylim(b, t) # update the ylim(bottom, top) values

    # ax = heatmap.ax_heatmap
    if is_imputed is not None:
        for i in range(is_imputed.shape[0]):
            if is_imputed[i,0]:
                for j in range(is_imputed.shape[1]):
                    # heatmap.add_patch(Rectangle((j, is_imputed.shape[0]-1-i), closed=True, fill=False, edgecolor='gray', lw=3))
                    box = np.array([[0, is_imputed.shape[0]-1-i], [4, is_imputed.shape[0]-1-i], [4, is_imputed.shape[0]-1-i+1], [0, is_imputed.shape[0]-1-i+1]])
                    heatmap.add_patch(Polygon(box, closed=True, fill=False, edgecolor='gray', lw=1.5, ls='--', clip_on=False))

    if output_path is not None:
        plt.savefig(output_path)
    else:
        plt.show()

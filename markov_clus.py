"""
This script creates gene families using markov clustering
"""
import pandas as pd
import numpy as np
import networkx as nx
import markov_clustering as mc
pd.options.display.width = 0

data = pd.read_csv('/nam-99/ablage/nam/peleke/result_TAIR10_filtered.csv')

# Generate graph with networkx
# query and subject genes will serve as  nodes and bitscores as the weight of edges
edges_with_weights = [(data['qgene'][i], data['sgene'][i], 
                       data['bitscore'][i]) for i in range(data.shape[0])]
G = nx.Graph()
G.add_weighted_edges_from(edges_with_weights)
nodes = list(G.nodes())

# Getting adjacency matrix
matrix = nx.to_scipy_sparse_matrix(G)

# MCL
results = mc.run_mcl(matrix, inflation=1.1)
hard_clusters = mc.get_clusters(results, keep_overlap=False)

# Produce families
gene_families = []
for i, tup in enumerate(hard_clusters):
    for idx in tup:
        gene_families.append([i, nodes[idx]])

# writing to csv
fam_df = pd.DataFrame(gene_families, columns=['family_id', 'Gene_id'])
fam_df.to_csv('/nam-99/ablage/nam/peleke/gene_families.csv', index=False)

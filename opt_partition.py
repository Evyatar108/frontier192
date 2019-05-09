from collections import defaultdict

import networkx as nx
import pickle

from sklearn import cluster


def __partition_by_labels(G, labels):
    biGraph = nx.Graph()
    node_to_label = {node: label for node, label in zip(G.nodes, labels)}

    clusters_dict = defaultdict(set)
    for node, label in node_to_label.items():
        clusters_dict[label].add(node)
    B = max(clusters_dict.values(), key=len)

    node_to_label = {node: 0 if node in B else 1 for node in G.nodes}

    for node, label in node_to_label.items():
        for adj_node in G.adj[node].keys():
            if node_to_label[adj_node] != node_to_label[node]:
                biGraph.add_edge(node, adj_node)

    B_bigraph = biGraph.subgraph(B)

    matching = nx.bipartite.maximum_matching(biGraph, top_nodes=B_bigraph)  # use hopcroft-karp algorithm

    nodes_to_remove = nx.bipartite.to_vertex_cover(biGraph, matching, top_nodes=B_bigraph)

    B = B.difference(nodes_to_remove)
    A = set(G.nodes).difference(B).difference(nodes_to_remove)

    return A, B, nodes_to_remove


def find_opt_partiton_kmeans(G):
    adj_mat = nx.to_scipy_sparse_matrix(G)
    kmeans = cluster.KMeans(n_clusters=2, n_init=10)
    kmeans.fit(adj_mat)
    return __partition_by_labels(G, kmeans.labels_)


def find_opt_partiton_spectral(G, clusters_num):
    adj_mat = nx.to_scipy_sparse_matrix(G)
    sc = cluster.SpectralClustering(clusters_num, affinity='precomputed', n_init=100)
    sc.fit(adj_mat)
    return __partition_by_labels(G, sc.labels_)


def find_opt_partiton_dbscan(G):
    adj_mat = nx.to_scipy_sparse_matrix(G)
    sc = cluster.DBSCAN(metric='precomputed')
    sc.fit(adj_mat)
    return __partition_by_labels(G, sc.labels_)


if __name__ == '__main__':
    pickle_in = open("nx_graph.pickle", "rb")
    G = pickle.load(pickle_in)
    find_opt_partiton(G)

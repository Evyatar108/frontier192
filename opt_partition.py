from collections import defaultdict
from itertools import chain
import kaHIP
import networkx as nx
import metis

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
    kmeans = cluster.KMeans(n_clusters=2, n_init=1)
    kmeans.fit(adj_mat)
    return __partition_by_labels(G, kmeans.labels_)


def find_opt_partition_spectral(G, clusters_num):
    adj_mat = nx.to_scipy_sparse_matrix(G)
    sc = cluster.SpectralClustering(clusters_num, affinity='precomputed', n_init=1, eigen_solver='amg')
    sc.fit(adj_mat)
    return __partition_by_labels(G, sc.labels_)

def find_opt_partition_cuts(G: nx.Graph):
    node_cuts = list(nx.all_node_cuts(G))
    nodes_to_remove = set(chain.from_iterable(node_cuts))
    V = set(G)
    AB_nodes = V.difference(nodes_to_remove)
    new_G = G.subgraph(AB_nodes)
    B_subgraph = max(nx.connected_component_subgraphs(new_G), key=len)
    B = set(B_subgraph)
    A = set(node for node in new_G if node not in B)
    return A, B, nodes_to_remove


def find_opt_partiton_dbscan(G):
    adj_mat = nx.to_scipy_sparse_matrix(G)
    sc = cluster.DBSCAN(metric='precomputed')
    sc.fit(adj_mat)
    return __partition_by_labels(G, sc.labels_)


def find_opt_partition_kahip(G: nx.Graph):
    node_list = list(G)
    node_to_index = {node_list[index]:index for index in range(len(node_list))}
    xadj = [0]
    sum = 0
    for v in node_list:
        n_neighbores = len(list(G.neighbors(v)))
        xadj.append(sum + n_neighbores)
        sum += n_neighbores

    adjncy = []
    for v in node_list:
        for u in G.neighbors(v):
            adjncy.append(node_to_index[u])

    num_of_deleted_nodes, deleted_nodes_indexes = kaHIP.node_separator(len(G), None, xadj, None, adjncy, 2, 0.03, True, 1, kaHIP.STRONG)
    deleted_nodes = set(node_list[deleted_node_index] for deleted_node_index in deleted_nodes_indexes)
    new_G = G.subgraph(set(G).difference(deleted_nodes))
    B_subgraph = max(nx.connected_component_subgraphs(new_G), key=len)
    B = set(B_subgraph)
    A = set(node for node in new_G if node not in B)
    return A, B, deleted_nodes

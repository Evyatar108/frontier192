import sys

import numpy as np

from opt_partition import find_opt_partiton_kmeans, find_opt_partiton_spectral, find_opt_partiton_dbscan

sys.path.insert(0, "/home/cluster/users/siditom/data/phd/KmerDB/aux/lib/python3.4/site-packages/networkx")
from approx_partition import find_approx_partition
import networkx as nx
import pickle

np.random.seed(1)
clusters_num = 10
print('clusters num each opt iteration: ' + clusters_num)

def save_object(V, path):
    pickle_out = open(path, "wb")
    pickle.dump(V, pickle_out)
    pickle_out.close()


def load_object(path):
    pickle_in = open(path, "rb")
    return pickle.load(pickle_in)


G = load_object("nx_graph.pickle")

def remove_edges_with_weight_of(G, mw):
    E = G.edges(data='weight')
    print("number of edges pre-filtering = " + str(len(E)))
    Ew = [(u, v, w) for (u, v, w) in E if w < mw]
    print("number of edges post-filtering = " + str(len(Ew)))
    G.remove_edges_from(Ew)
    return G


def graph_vsize_attd(G):
    V = set(G)
    for vi in V:
        G.nodes[vi]['size'] = 1
    return G


def remove_double_edges(G):
    newG = nx.Graph()
    newG.add_nodes_from(set(G))
    newG.add_edges_from(set(G.edges()))
    for (u, v) in newG.edges():
        newG[u][v]['weight'] = 1
    return newG


def partition(G):
    A = set()
    B_subgraph = G
    B = set(B_subgraph)
    new_B = B
    deleted_nodes = set()
    thresh = 0.5
    ratio = 0.1
    while ratio < thresh:
        deleted_nodes_ext = G.nodes()
        A_ext = set()
        for i in range(1):
            temp_A_ext, temp_B, temp_deleted_nodes_ext = find_opt_partiton_spectral(B_subgraph)
            if len(deleted_nodes_ext) > len(temp_deleted_nodes_ext):
                deleted_nodes_ext = temp_deleted_nodes_ext
                A_ext = temp_A_ext
                new_B = temp_B
        B = new_B
        B_subgraph = G.subgraph(B)
        B_connected_components = list(nx.connected_components(B_subgraph))
        print('B has '+str(len(B_connected_components))+' connected components, moving the small ones to A')
        B = max(B_connected_components, key=len)
        B_connected_components.remove(B)
        B_subgraph = G.subgraph(B)
        A_ext_from_connected_components = set()
        for B_connected_component in B_connected_components:
            A_ext_from_connected_components.update(B_connected_component)

        A.update(A_ext)
        A.update(A_ext_from_connected_components)
        deleted_nodes.update(deleted_nodes_ext)

        ratio = len(A) / (len(B) + len(A))
        print('current deleted nodes: '+ str(len(deleted_nodes)))
        print('current deleted nodes ext: '+ str(len(deleted_nodes_ext)))
        print('Current |A_ext| = ' + str(len(A_ext)))
        print('Current |A_ext_from_connected_components| = ' + str(len(A_ext_from_connected_components)))
        print('Current |A| = ' + str(len(A)))
        print('Current |B| = ' + str(len(B)))
        print('ratio = ' + str(ratio))
    return A, B, deleted_nodes

G = remove_edges_with_weight_of(G, 0.3)
G = G.to_undirected()
G = remove_double_edges(G)
G = graph_vsize_attd(G)
connected_components_subgraphs = nx.connected_component_subgraphs(G)
print("Number of |G.V| = " + str(len(G)) + ", Number of |G.E| = " + str(len(G.edges)))
nV = len(G)

print("Calculated connected components")
nodes_clusters = []
deleted_nodes = []
for connected_component_subgraph in connected_components_subgraphs:
    if len(connected_component_subgraph) > (nV / 2):
        A, B, deleted_nodes = partition(connected_component_subgraph)
        print("END Cut:")
        print("Group |A| = " + str(len(A)))
        print("Group |B| = " + str(len(B)))
        print("Group |Deleted| = " + str(len(deleted_nodes)))
        nodes_clusters.append(A)
        nodes_clusters.append(B)
    else:
        nodes_clusters.append(connected_component_subgraph.nodes())

save_object(deleted_nodes, "deleted_nodes.pickle")

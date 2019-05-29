import numpy as np
from opt_partition import find_opt_partition_kahip, find_opt_partition_spectral
import networkx as nx
import pickle
from pyamg import smoothed_aggregation_solver

clusters_num = 2
print('clusters num each opt iteration: ' + str(clusters_num))


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


def try_load_cached_partition(path):
    try:
        A_cached, B_cached, deleted_nodes_cached = load_object(path)
    except:
        A_cached, B_cached, deleted_nodes_cached = (set(), set(G), set(G))
        save_object((A_cached, B_cached, deleted_nodes_cached), path)

    return A_cached, B_cached, deleted_nodes_cached


def partition(G):
    thresh = 0.5
    ratio = 0.1

    A, B, deleted_nodes = find_opt_partition_kahip(G)
    A_cached, B_cached, deleted_nodes_cached = try_load_cached_partition('kahip_cache.pickle')
    if len(deleted_nodes) < len(deleted_nodes_cached):
        save_object((A, B, deleted_nodes), 'kahip_cache.pickle')
        print('Improved cached kahip from ' + str(len(deleted_nodes_cached)) + ' deleted nodes to ' + str(
            len(deleted_nodes)) + ' deleted nodes')
    else:
        A, B, deleted_nodes = A_cached, B_cached, deleted_nodes_cached
        print('Loaded cached kahip partition, no improvement found')

    ratio = len(A) / (len(B) + len(A))
    new_B = B
    B_subgraph = G.subgraph(B)

    print('Current deleted nodes: ' + str(len(deleted_nodes)))
    print('Current |A| = ' + str(len(A)))
    print('Current |B| = ' + str(len(B)))
    print('ratio = ' + str(ratio))

    while False: # ratio < thresh:
        deleted_nodes_ext = G.nodes()
        A_ext = set()
        for i in range(10):
            temp_A_ext, temp_B, temp_deleted_nodes_ext = find_opt_partition_spectral(B_subgraph, clusters_num)
            if len(temp_A_ext) / len(temp_deleted_nodes_ext) > len(A_ext) / len(deleted_nodes_ext):
                deleted_nodes_ext = temp_deleted_nodes_ext
                A_ext = temp_A_ext
                new_B = temp_B
        B = new_B
        B_subgraph = G.subgraph(B)
        B_connected_components = list(nx.connected_components(B_subgraph))
        print('B has ' + str(len(B_connected_components)) + ' connected components, moving the small ones to A')
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
        print('Current deleted nodes: ' + str(len(deleted_nodes)))
        print('Current deleted nodes ext: ' + str(len(deleted_nodes_ext)))
        print('Current |A_ext| = ' + str(len(A_ext)))
        print('Current |A_ext_from_connected_components| = ' + str(len(A_ext_from_connected_components)))
        print('Current |A| = ' + str(len(A)))
        print('Current |B| = ' + str(len(B)))
        print('ratio = ' + str(ratio))

    A_cached, B_cached, deleted_nodes_cached = try_load_cached_partition('end_cache.pickle')
    if len(deleted_nodes) < len(deleted_nodes_cached):
        save_object((A, B, deleted_nodes), 'end_cache.pickle')
        print('Improved cached end partition from ' + str(len(deleted_nodes_cached)) + ' deleted nodes to ' + str(
            len(deleted_nodes)) + ' deleted nodes')
    else:
        print('End partition loaded, no improvement found')
        A, B, deleted_nodes = A_cached, B_cached, deleted_nodes_cached

    return A, B, deleted_nodes


G = remove_edges_with_weight_of(G, 0.3)
G = G.to_undirected()
G = remove_double_edges(G)
G = graph_vsize_attd(G)
connected_components_subgraphs = nx.connected_component_subgraphs(G)
print("Number of |G.V| = " + str(len(set(G))) + ", Number of |G.E| = " + str(len(set(G.edges))))
nodes_clusters = []

while True:
    A, B, deleted_nodes = partition(G)
    print("END Cut:")
    print("Group |A| = " + str(len(A)))
    print("Group |B| = " + str(len(B)))
    print("Group |Deleted| = " + str(len(deleted_nodes)))

import numpy as np
# sys.path.insert(0, "/home/cluster/users/siditom/data/phd/KmerDB/aux/lib/python3.4/site-packages/networkx")
from approx_partition import partition_approx
import networkx as nx
import pickle

from opt_partition import find_opt_partiton

np.random.seed(1)

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
    Ew = [(u, v, w) for (u, v, w) in E if w < mw]
    print("number of Edges = " + str(len(Ew)))
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


G = remove_edges_with_weight_of(G, 0.3)
G = G.to_undirected()
G = remove_double_edges(G)
G = graph_vsize_attd(G)
connected_components_subgraphs = nx.connected_component_subgraphs(G)
print("Number of |G.V| = " + str(len(G)) + ", Number of |G.E| = " + str(len(G.edges)))
nV = len(G)

print("Calculated connected components")
nodes_clusters = []
for connected_component_subgraph in connected_components_subgraphs:
    if len(connected_component_subgraph) > (nV / 2):
        A, B, deleted_nodes = find_opt_partiton(connected_component_subgraph)
        print("END Cut:")
        print("Group |A| = " + str(len(A)))
        print("Group |B| = " + str(len(B)))
        print("Group |Deleted| = " + str(len(deleted_nodes)))
        nodes_clusters.append(A)
        nodes_clusters.append(B)
    else:
        nodes_clusters.append(connected_component_subgraph.nodes())

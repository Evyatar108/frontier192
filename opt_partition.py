import networkx as nx
from sklearn.cluster import SpectralClustering
import pickle


def find_opt_partiton(G: nx.DiGraph):
    adj_mat = nx.to_numpy_matrix(G)

    sc = SpectralClustering(2, affinity='precomputed', n_init=100)
    sc.fit(adj_mat)

    print('spectral clustering')

    biGraph = nx.Graph()
    node_to_label = {node: label for node, label in zip(G.nodes, sc.labels_)}

    for node, label in node_to_label.items():
        for adj_node in G.adj[node].keys():
            if node_to_label[adj_node] != node_to_label[node]:
                biGraph.add_edge(node, adj_node)

    matching = nx.bipartite.maximum_matching(biGraph) # use hopcroft-karp algorithm

    nodes_to_remove = nx.bipartite.to_vertex_cover(biGraph, matching)

    G.remove_nodes_from(nodes_to_remove)

    A = [node for node in G.nodes if node_to_label[node] == 0]
    B = [node for node in G.nodes if node_to_label[node] == 1]

    return A, B, nodes_to_remove


if __name__ == '__main__':
    pickle_in = open("nx_graph.pickle", "rb")
    G = pickle.load(pickle_in)
    find_opt_partiton(G)

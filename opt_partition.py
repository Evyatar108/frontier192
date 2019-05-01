import numpy as np
import networkx as nx
from sklearn.cluster import SpectralClustering
import pickle
np.random.seed(1)


def find_opt_partiton(G: nx.DiGraph):
    adj_mat = nx.to_numpy_matrix(G)

    sc = SpectralClustering(2, affinity='precomputed', n_init=100)
    sc.fit(adj_mat)

    print('spectral clustering')
    # print(sc.labels_)
    # print(f'initial graph nodes: {test_graph.G.nodes}')

    biGraph = nx.Graph()
    node_to_label = {node: label for node, label in zip(G.nodes, sc.labels_)}

    for node, label in node_to_label.items():
        for adj_node in G.adj[node].keys():
            if node_to_label[adj_node] != node_to_label[node]:
                biGraph.add_edge(node, adj_node)

    # print(f'bipartite graph created from clustering: {biGraph.nodes}')

    matching = nx.bipartite.maximum_matching(biGraph) # use hopcroft-karp algorithm
    # print(f'bipartite graph maximum matching: {matching}')

    nodes_to_remove = nx.bipartite.to_vertex_cover(biGraph, matching)
    # print(f'nodes to remove: {nodes_to_remove}')

    G.remove_nodes_from(nodes_to_remove)
    # print(f'new graph (after deletion): {G.nodes}')

    left = [node for node in G.nodes if node_to_label[node] == 0]
    right = [node for node in G.nodes if node_to_label[node] == 1]
    print(f'size of left group: {len(left)}')
    print(f'size of right group: {len(right)}')
    print(f'number of deleted nodes: {len(nodes_to_remove)}')


if __name__ == '__main__':
    pickle_in = open("nx_graph.pickle", "rb")
    G = pickle.load(pickle_in)
    find_opt_partiton(G)

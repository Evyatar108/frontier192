#sys.path.insert(0, "/home/cluster/users/siditom/data/phd/KmerDB/aux/lib/python3.4/site-packages/networkx")
import networkx as nx
import numpy as np

def __remove_cut_min_vertex_cover(G, A, B):

    biGraph = nx.Graph()

    newB = set()
    for vi in A:
        newB.update(B.intersection(G.neighbors(vi)))
    node_to_label_B = {node: 0 for node in newB}

    newA = set()
    for vi in B:
        newA.update(A.intersection(G.neighbors(vi)))
    node_to_label_A = {node: 1 for node in newA}

    node_to_label_A.update(node_to_label_B)

    for node, label in node_to_label_A.items():
        for adj_node in G.adj[node].keys():
            if adj_node in node_to_label_A and node_to_label_A[adj_node] != node_to_label_A[node]:
                biGraph.add_edge(node, adj_node)

    matching = nx.bipartite.maximum_matching(biGraph, newA)  # use hopcroft-karp algorithm

    nodes_to_remove = nx.bipartite.to_vertex_cover(biGraph, matching, newA)

    A = A.difference(nodes_to_remove)
    B = B.difference(nodes_to_remove)

    print(f'size of newA group: {len(A)}')
    print(f'size of newB group: {len(B)}')
    print(f'number of deleted nodes: {len(nodes_to_remove)}')

    return A, B, nodes_to_remove


def find_approx_partition(G):
    print("Partition connected component")
    V = set(G)
    A = set(np.random.choice(G, int((len(V) / 3) * 2)))
    B = V.difference(A)
    print("Random cut separation")
    print("Iteration to modify cut:")

    transferred_nodes = 1
    iteration = 0
    while transferred_nodes >= 1:
        transferred_nodes = 0
        print("Iteration = " + str(iteration))
        iteration += 1

        print("Group |A| = " + str(len(A)))
        print("Group |B| = " + str(len(B)))

        for vi in A.copy():
            all_neighbors = set(G.neighbors(vi))
            A_neighbors_num = len(all_neighbors.difference(B))
            B_neighbors_num = len(all_neighbors.difference(A))
            if A_neighbors_num < B_neighbors_num:
                A.remove(vi)
                B.add(vi)
                transferred_nodes += 1

        for vi in B.copy():
            all_neighbors = set(G.neighbors(vi))
            A_neighbors_num = len(all_neighbors.difference(B))
            B_neighbors_num = len(all_neighbors.difference(A))
            if B_neighbors_num < A_neighbors_num:
                B.remove(vi)
                A.add(vi)
                transferred_nodes += 1

    print("Finished modifying cut")
    print("Cut Size:")
    print("New group |A| = " + str(len(A)))
    print("New group |B| = " + str(len(B)))

    print("Remove neighboring nodes")

    A, B, deleted_nodes = __remove_cut_min_vertex_cover(G, A, B)
    return A, B, deleted_nodes
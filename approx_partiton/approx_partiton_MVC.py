import networkx as nx


def remove_cut_min_vertex_cover(G, A, B):

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

    # print(f'the bipartite graph that created from clustering: {biGraph}')

    matching = nx.bipartite.maximum_matching(biGraph, newA)  # use hopcroft-karp algorithm

    nodes_to_remove = nx.bipartite.to_vertex_cover(biGraph, matching, newA)

    G.remove_nodes_from(nodes_to_remove)
    A = A.difference(nodes_to_remove)
    B = B.difference(nodes_to_remove)

    print(f'size of newA group: {len(A)}')
    print(f'size of newB group: {len(B)}')
    print(f'number of deleted nodes: {len(nodes_to_remove)}')

    return A, B, nodes_to_remove
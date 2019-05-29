from networkx import DiGraph

from opt_partition import find_opt_partition_kahip, find_opt_partition_spectral

G = DiGraph()
G.add_edge(6, 2, weight=0.51)
G.add_edge(6, 8, weight=0.45)
G.add_edge(3, 8, weight=0.49)
G.add_edge(5, 1, weight=0.4679)
G.add_edge(2, 5, weight=0.54858)
G.add_edge(1, 4, weight=0.467658)
G.add_edge(1, 2, weight=0.56776)
G.add_edge(2, 3, weight=0.667)
G.add_edge(4, 3, weight=0.567568)
G.add_edge(3, 7, weight=0.487884)
G.add_edge(5, 4, weight=0.56767)
G.add_edge(7, 9, weight=0.65567)
G.add_edge(10, 7, weight=0.576367)
G.add_edge(10, 6, weight=0.7878)
G.add_edge(14, 9, weight=0.64887)
G.add_edge(10, 17, weight=0.56756)
G.add_edge(10, 15, weight=0.67876)
G.add_edge(10, 14, weight=0.468747)
G.add_edge(15, 14, weight=0.56685)
G.add_edge(11, 16, weight=0.6577)
G.add_edge(17, 12, weight=0.56776)
G.add_edge(13, 17, weight=0.57676)
G.add_edge(14, 11, weight=0.6758)
G.add_edge(12, 13, weight=0.6585)
G.add_edge(12, 11, weight=0.5887)
G.add_edge(7, 4, weight=0.5678)
G.add_edge(8, 7, weight=0.67868)
G.add_edge(17, 9, weight=0.5986)

nodes_to_remove = find_opt_partition_kahip(G.to_undirected())
print(f'nodes to remove = {nodes_to_remove}')

A, B, nodes_to_remove = find_opt_partition_spectral(G.to_undirected(), 2)
print(f'nodes to remove = {nodes_to_remove}')
print('|A| = ' + str(len(A)) + ', |B| = ' + str(len(B)) + ', |nodes_to_remove| = ' + str(len(nodes_to_remove)))






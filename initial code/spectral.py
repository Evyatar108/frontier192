import numpy as np
import networkx as nx
from sklearn.cluster import SpectralClustering
from sklearn import metrics
from graph_constraction import load_object
np.random.seed(1)


G = load_object("nx_graph.pickle")
adj_mat = nx.to_numpy_matrix(G)

sc = SpectralClustering(2, affinity='precomputed', n_init=100)
sc.fit(adj_mat)


print('spectral clustering')
print(sc.labels_)

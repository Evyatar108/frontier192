import numpy as np
import networkx as nx
from sklearn.cluster import SpectralClustering
from sklearn import metrics
import pickle
np.random.seed(1)


pickle_in = open("nx_graph.pickle", "rb")
G = pickle.load(pickle_in)
adj_mat = nx.to_numpy_matrix(G)

sc = SpectralClustering(2, affinity='precomputed', n_init=100)
sc.fit(adj_mat)


print('spectral clustering')
print(sc.labels_)

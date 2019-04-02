import sys
import os
import time
import numpy as np
import datetime
import shutil
from subprocess import Popen
from subprocess import call
sys.path.insert(0, "/home/cluster/users/siditom/data/phd/KmerDB/aux/lib/python3.4/site-packages/networkx")
import networkx as nx
import networkx.utils 
import glob
import re, string
import pickle
import random 
from networkx.algorithms.flow import shortest_augmenting_path
from sklearn.cluster import SpectralClustering
from sklearn import metrics
import operator
V = []

def get_vertices():
	path = "/home/cluster/users/siditom/data/phd/KmerDB/data/pdb_chains/blast_allVSall_Chain/nr100/nr100.fasta"
	with open(path,"r") as f:
		v_all = [line[1:] for line in f.read().split("\n") if line.startswith(">")]
	V = dict()
	for vi in v_all:
		V[vi] = list()
	return V

def get_edges(V):
	path = "/home/cluster/users/siditom/data/phd/KmerDB/data/pdb_chains/blast_allVSall_Chain/"
	avsa_blast = [f for f in glob.glob(path+"*.out")]
	
	inx=1
	s = len(avsa_blast)
	for f in avsa_blast:
		with open(f,"r") as ff:
			edges_text = [l.split("\t") for l in ff.readlines()]
			edges = [[e[0].replace(">",""),e[1].replace("_",""),{"weight": float(e[2])}] for e in edges_text]
		#	edges = [[e[0].replace(">",""),e[1].replace("_",""), int(float(e[2]))] for e in edges_text]
		for e in edges:
			if (e[0] in V) and (e[1] in V) and (e[0] != e[1]):
				V[e[0]].append(e)
#		print ("SplitPdbData: added edges. blast file " + str(inx) + "/" + str(s))
		inx = inx + 1
	
	return V

def save_object(V,path):

	pickle_out = open(path,"wb")
	pickle.dump(V, pickle_out)
	pickle_out.close()


def load_object(path):
	pickle_in = open(path,"rb")
	return pickle.load(pickle_in)	

def constract_graph(V):
	G = nx.DiGraph()
#	G = nx.Graph()
	G.add_nodes_from(V.keys())
	inx = 0
	s = len(V.keys())
	for vi in V:
		#es = [e for e in V[vi]]
		es = [[e[0],e[1],(e[2]["weight"]/100)] for e in V[vi]]
		#G.add_edges_from(es)			
		G.add_weighted_edges_from(es)			
		print ("SplitPdbData: added edges. blast file " + str(inx) + "/" + str(s))
		inx = inx+1
	return G

def generate_net(G):
	V = list(G)
	G.add_node("s")
	G.add_node("t")
	A = np.random.choice(V,10000)
	for vi in A:
		#G.add_weighted_edges_from([["s",vi,float(200)]])
		G.add_weighted_edges_from([[vi,"t",float(100)]])
	A = np.random.choice(V,100)
	for vi in A:
		G.add_weighted_edges_from([["s",vi,float(100)]])
	#E = G.edges(data = 'weight')
	E = G.edges()
	for (u,v) in E:
		if (u != "s") and (v != "t"):
			G[u][v]['capacity'] = 1
		else:
			G[u][v]['capacity'] = 100
			
		if (not (v,u) in E) and (u != "s") and (v != "t"):
			G.add_weighted_edges_from([[v,u,1]])
			G[v][u]['capacity'] = 1
	return G



def graph_partition(G):

	#cut_value, partition = nx.minimum_cut(G,"s","t",capacity='weight',flow_func=shortest_augmenting_path)
	cut_value, partition = nx.minimum_cut(G,"s","t",flow_func=shortest_augmenting_path)

	reachable, non_reachable = partition
	print(len(reachable))
	print(len(non_reachable))
	print(cut_value)

######################################
######### Get Raw Graph Data #########
######################################
"""
V = get_vertices()
print("Recovered Vertices ... ")
V = get_edges(V)
print("Recovered Edges ... ")
save_object(V,"raw_graph.pickle")
"""

######################################
######### Load Raw Graph Data ########
######################################
"""
V = load_object("raw_graph.pickle")
"""


######################################
############ Create Graph ############
######################################
"""
G = constract_graph(V)
save_object(G,"nx_graph.pickle")
"""

######################################
############# Load Graph #############
######################################
G = load_object("nx_graph.pickle")


######################################
############# Calculate  #############
######################################

with open("edges.dat","w") as f:
#	f.write(str(G.degree).replace(",","\n"))
	f.write(str( sorted(G.edges(data = 'weight'),key=lambda x:x[2],reverse=True)))
with open("in.dat","w") as f:
#	f.write(str(G.degree).replace(",","\n"))
	f.write(str( sorted(G.in_degree,key=lambda x:x[1],reverse=True)))
with open("out.dat","w") as f:
#	f.write(str(G.degree).replace(",","\n"))
	f.write(str(G.out_degree))


"""
H = nx.Graph(G)
IS = nx.maximal_independent_set(H)
print(len(IS))


ccs = nx.weakly_connected_components(G)
"""

def remove_edges_with_weight_of(G,mw):
	E = G.edges(data='weight')
	Ew = [(u,v,w) for (u,v,w) in E if w < mw]
	print("number of Edges = "+str(len(Ew)))
	G.remove_edges_from(Ew)
	return G

def graph_vsize_attd(G):
	V = set(G)
	for vi in V:
		G.nodes[vi]['size']=1
	return G

def combine_group_in_graph(G, s):
	A = set(G.neighbors(s))
	c_size = G.nodes[s]['size']

	for vi in A:
		c_size+=G.nodes[vi]['size']
	G.nodes[s]['size'] = c_size
	#print("Node "+str(s)+" , size = "+str(G.nodes[s]['size']))
	if not A.issubset(set(G)):
		raise Exception
	A_neighbors = set()
	for v in A:
		A_neighbors.update(set(G.neighbors(v)))
	A_neighbors.discard(s)
	A_neighbors = A_neighbors.difference(A)
	G.remove_nodes_from(A)
	sE = [[s,v] for v in A_neighbors]
	G.add_edges_from(sE)
	return G

def spectral(c):
	print("Spectral Clustering ...")
	adj_mat = nx.to_numpy_matrix(c)
	sc = SpectralClustering(2, affinity='precomputed', n_init=100)
	sc.fit(adj_mat)
	print('spectral clustering')
	print(sc.labels_)
	print(np.sum(sc.labels_))
	return sc.labels_

def remove_double_edges(G):
	newG = nx.Graph()
	V = set(G)
	newG.add_nodes_from(V)
	E = set(G.edges())
	newG.add_edges_from(E)
	newE = newG.edges()
	
	for (u,v) in newE:
		newG[u][v]['weight'] = 1
	return newG

def partition(G):
	print("Partition Graph")
	G = remove_double_edges(G)
	V = set(G)
	print("|V| = "+str(len(V)))
	A = set(np.random.choice(G,int((len(V)/3) * 2)))
	B = V.difference(A)
	print("Random cut separation")
	#print("Cut size = "+str(nx.cut_size(G, A, B,weight='weight')))
	#Iterations to change the random cut:
	print("Iteration to modify cut:")
	for i in range(1,10):
		print("Iteration = "+str(i))
		ABneighbors = set()
		for vi in A:
			all_neighbors = set(G.neighbors(vi))
			B_neighbors = all_neighbors.difference(A)
			ABneighbors.update(B_neighbors)
		BAneighbors = set()
		
		for vi in B:
			all_neighbors = set(G.neighbors(vi))
			A_neighbors = all_neighbors.difference(B)
			BAneighbors.update(A_neighbors)
		print("Cut size A to B = "+str(len(ABneighbors)))
		print("Cut size B to A = "+str(len(BAneighbors)))
		print("Group |A| = "+str(len(A)))
		print("Group |B| = "+str(len(B)))

		#currA = set(np.random.choice(list(A),100))
		#currB = set(np.random.choice(list(B),100))
		currA = A.copy()
		currB = B.copy()
		scoreA = dict()
		scoreB = dict()
		for vi in currA:
			all_neighbors = set(G.neighbors(vi))
			A_neighbors_num = len(all_neighbors.difference(B))
			B_neighbors_num = len(all_neighbors.difference(A))
			if (A_neighbors_num < B_neighbors_num):
				#scoreA[vi] = A_neighbors_num/B_neighbors_num 
				A.remove(vi)
				B.add(vi)
		for vi in currB:
			all_neighbors = set(G.neighbors(vi))
			A_neighbors_num = len(all_neighbors.difference(B))
			B_neighbors_num = len(all_neighbors.difference(A))
			if (B_neighbors_num < A_neighbors_num):
				#scoreB[vi] = B_neighbors_num/A_neighbors_num 
				B.remove(vi)
				A.add(vi)
		"""
		va = max(scoreA.items(), key=operator.itemgetter(1))[0]	
		A.remove(va)
		B.add(va)
		vb = max(scoreB.items(), key=operator.itemgetter(1))[0]	
		B.remove(vb)
		A.add(vb)
		"""
	print("Finished modifying cut")
	print("Cut Size:")
	print("Group |A| = "+str(len(A)))
	print("Group |B| = "+str(len(B)))

	print("Remove neighboring nodes")
	newB = B.copy()
	for vi in A:
		newB = newB.difference(G.neighbors(vi))
	dN = B.difference(newB)
	print("Post Remove Cut:")
	print("Group |A| = "+str(len(A)))
	print("Group |B| = "+str(len(B)))


	return A,newB,dN
	"""
	iteration=1
	while nx.cut_size(G, A, B,weight='weight') > 0 :
		print("Iteration = "+str(iteration))
		score = dict()
		neighbors = set()
		for vi in A:
			all_neighbors = set(G.neighbors(vi))
			B_neighbors = all_neighbors.difference(A)
			neighbors.update(B_neighbors)
			B_neighbors_num = len(B_neighbors)
			score[vi] = B_neighbors_num 
		for vi in B:
			all_neighbors = set(G.neighbors(vi))
			A_neighbors_num = len(all_neighbors.difference(B))
			score[vi] = A_neighbors_num 
		vi = max(score.items(), key=operator.itemgetter(1))[0]	
		A.discard(vi)
		B.discard(vi)
		print("Cut size = "+str(len(neighbors)))
		iteration += 1
	"""	


def calc_dict(c,A):
	c1 = c.subgraph(A)
			
	ccsc = [set(g) for g in nx.connected_component_subgraphs(c1)]
	vDict = dict()	
	for cg in ccsc:
		for vi in cg:
			vDict[vi] = len(cg)
	return c1,vDict

def regroup(c,A,dN):
	#B = set(np.random.choice(list(dN),i*1000))
	#A.update(B)

	c1,vDict = calc_dict(c,A)
	print(c1)
	print(vDict)
	for vi in dN:
		print(vi)
		neighbors_vi = set(c.neighbors(vi)).intersection(set(c1))
		s = 1
		if (neighbors_vi == None):
			continue
		for vj in neighbors_vi:
			s+=vDict[vj]
		if (s < len(A)/10):
			print(s)
			A.add(vi)
			c1,vDict = calc_dict(c,A)
			
	
	print("Size of wccs 1:")
	print(str([len(g) for g in ccsc]))
	return A


	
#G = generate_net(G)
#graph_partition(G)

		

G = remove_edges_with_weight_of(G,0.3)
ccs = nx.weakly_connected_component_subgraphs(G)
print("Number of nodes = "+str(len(G)) + ", Number of Edges = "+str(len(G.edges)))
nV = len(G)

nodes = set()
print("Calculated Weak Components")
for c in ccs:
	if len(c)>(nV/2):
	#	S = np.random.choice(list(c), int(len(c)/2) )	
			
		c = c.to_undirected()
		c = remove_double_edges(c)
		c = graph_vsize_attd(c)
		spectral(c)
		"""	
		print("Transformed to undirected graph")
		in_d = sorted(c.degree,key=lambda x:x[1],reverse=True)
		print("Calculated in degrees of all nodes")
		B = [v for (v,deg) in in_d]
		print("Number of nodes = "+str(len(c)))
		B = B[1:1000]
		for vi in B:
			if (vi in c):
				c = combine_group_in_graph(c,vi)
				print("combine_grop_in_graph: c.nodes["+vi+"]['size'] = "+str(c.nodes[vi]['size']))
				#print("Number of nodes = "+str(len(c)) + ", Number of Edges = "+str(len(c.edges)))
		print("Number of nodes = "+str(len(c)))
		"""

		"""
		A = set()
		dN = set()
		thresh = 0.5
		ratio = 0.1
		c1 = c.copy()
		while ratio < thresh :
			D,B,deleted_N = partition(c1)
			A.update(D)
			dN.update(deleted_N)
			c1 = c1.subgraph(B)
		
			print("Group |A| = "+str(len(A)))
			print("Group |B| = "+str(len(B)))
			print("Group |Deleted| = "+str(len(dN)))
			
			ratio = len(A)/(len(B)+len(A))
		"""


		"""	
		print("END Cut:")
		print("Group |A| = "+str(len(A)))
		print("Group |B| = "+str(len(B)))
		print("Group |Deleted| = "+str(len(dN)))
		A.update(B)
		c1 = c.subgraph(A)
			
		ccsc = nx.connected_component_subgraphs(c1)
		
		print("Size of wccs 0:")
		print(str([len(g) for g in ccsc]))
		
		A = regroup(c,A,dN)
		"""


		"""
		print("Actual size (including combined nodes):")
		lsize = list()
		for gc in ccsc:
			l = 0
			for v in gc:
				l += c.nodes[v]['size']
			lsize.append([len(gc),l])
		print(lsize)
		"""



"""		
		print("####")
		save_object(c,"reduced_1d_cluster.pickle")	
		with open("Edges20_.csv","w") as f:
			E = c.edges()
			for (u,v) in E:
				f.write(str(u)+","+str(v)+"\n")	

		#c = generate_net(c)
		#print("Net Generated")
		#graph_partition(c)
	
		nodes.update(A)
	else:
		nodes.update(set(c))


ccs = nx.weakly_connected_component_subgraphs(G.subgraph(nodes))
with open("wccg.dat","w") as f:
	for c in ccs:
		s = ""
		for vi in c:
			s += str(vi) + " "
		f.write(s+"\n")
	
"""
import numpy as np
import networkx as nx
from karateclub import Graph2Vec
from karateclub import GL2Vec
import os
import pandas as pd

def py_get_graph2Vec_embedding(graphsPath,wl_it=4,ebd_dim=8,attrib=False,base_feat=False,down_samp=0.0001,n_epochs=100,lr=0.015):	
	graphs = []
	files= os.listdir(graphsPath)
	for file in files:
		
		G = nx.read_graphml(graphsPath+file)
		G = nx.convert_node_labels_to_integers(G)
		graphs.append(G)

	print('prepared graphs list')

	graph2Vec_model = Graph2Vec(wl_iterations=wl_it, attributed=attrib, erase_base_features = not base_feat, dimensions=ebd_dim, workers=4,
					down_sampling=down_samp, epochs=n_epochs, learning_rate=lr)

	graph2Vec_model.fit(graphs)

	print('Graph2Vec fitted')

	embeddings = graph2Vec_model.get_embedding()
	
	names = []
	for fi in files:
		names.append(fi.split('.')[0])
	df = pd.DataFrame(data=embeddings)
	df['name'] = names
	
	return df
	
def py_get_graph2Vec_rectoverso_embedding(graphsPath,wl_it=4,ebd_dim=8,attrib=False,base_feat=False,down_samp=0.0001,n_epochs=100,lr=0.015):	
	graphs = []
	files= os.listdir(graphsPath)
	for file in files:
		G = nx.read_graphml(graphsPath+file)
		G = nx.convert_node_labels_to_integers(G)
		graphs.append(G)

	print('prepared graphs list')
	
	# Make a Graph2Vec pass for the initial edges direction 
	graph2Vec_model = Graph2Vec(wl_iterations=wl_it, attributed=attrib, erase_base_features = not base_feat, dimensions=ebd_dim, workers=4,
					down_sampling=down_samp, epochs=n_epochs, learning_rate=lr)

	graph2Vec_model.fit(graphs)
	print('Graph2Vec fitted for initial edges direction')

	embeddings = graph2Vec_model.get_embedding()

	# Reverse each directed graph
	i=0
	for G in graphs:
		G_ = G.reverse(copy=True)
		graphs[i] = G_
		i=i+1

	# Make another pass for the reverse edges direction 
	graph2Vec_model = Graph2Vec(wl_iterations=wl_it, attributed=attrib, dimensions=ebd_dim, workers=4,
					down_sampling=down_samp, epochs=n_epochs, learning_rate=lr)

	graph2Vec_model.fit(graphs)
	print('Graph2Vec fitted for reverse edges direction')

	embeddings2 = graph2Vec_model.get_embedding()
	embeddings = np.column_stack((embeddings,embeddings2))

	names = []
	for fi in files:
		names.append(fi.split('.')[0])

	df = pd.DataFrame(data=embeddings)
	df['name'] = names
	
	return df


def py_get_GL2Vec_embedding(graphsPath,wl_it=4,ebd_dim=8,attrib=False,down_samp=0.0001,n_epochs=100,lr=0.015):	
	graphs = []
	files= os.listdir(graphsPath)
	for file in files:
		G = nx.read_graphml(graphsPath+file)
		G = nx.convert_node_labels_to_integers(G)
		graphs.append(G)
	print('prepared graphs list')
	
	model = GL2Vec(wl_iterations=wl_it, dimensions=ebd_dim , workers=4,
				down_sampling=down_samp, epochs=n_epochs, learning_rate=lr)
	model.fit(graphs)
	print('GL2Vec fitted')

	embeddings = model.get_embedding()

	names = []
	for fi in files:
		names.append(fi.split('.')[0])
	df = pd.DataFrame(data=embeddings)
	df['name'] = names
	
	return df


#graphsDir = 'C:/Users/admbotella/Documents/pCloud local/boulot/data/Interaction Networks/tetrapods alpes/Graphs/'
#ebd = py_get_graph2Vec_embedding(graphsPath=graphsDir,wl_it=6,ebd_dim=32,attrib=False,down_samp=0.0001,n_epochs=100,lr=0.015)



import time
import os
from nltk import edit_distance
from pyclustering.cluster.kmedoids import kmedoids
import numpy as np
import pandas as pd
import argparse
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import networkx as nx
import matplotlib.pyplot as plt


import zipfile
with zipfile.ZipFile("diamond_distance_matrix.zip","r") as zip_ref:
    zip_ref.extractall(".")

mat = np.load('diamond_distance_matrix.npy')
mat = mat + mat.T

## get the adjaceny matrix ###
adj = np.zeros(mat.shape)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i,j]==1e+100:
            adj[i,j]=0
        else:
            adj[i,j]=1
        if i==j:
            adj[i,j]=0

## get the connected components ##
G = nx.from_numpy_matrix(adj)
sort=sorted(nx.connected_components(G), key = len, reverse=True)

## get AMP prototypes ###
prototype_id_all = []
for i in range(len(sort)):
    temp_mat = mat[list(sort[i]), :]
    temp_mat = temp_mat[:, list(sort[i])]
    ## get medoid for this connected components
    init_medoids = np.random.randint(temp_mat.shape[0], size=1)
    km = kmedoids(temp_mat, init_medoids, data_type='distance_matrix')
    km.process()
    centers = km.get_medoids()
    # get the id of center wrt the distance matrix
    prototype_id = list(sort[i])[centers[0]]
    prototype_id_all.append(prototype_id)

### get the corresponding AMPs from the prototype id ###
dataset = pd.read_csv('D1.csv')
dataset.keys()
amp = dataset.loc[dataset['Label'] == 1]
amp = amp.Sequences.tolist()

text_file = open("amp_prototype_diamond_distance.txt", "w")
for i in prototype_id_all:
    amp_prototype = amp[i]
    text_file.write(str(amp_prototype) + '\n')
text_file.close()

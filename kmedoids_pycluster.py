import numpy as np
from copy import deepcopy
import pandas as pd
from nltk.metrics.distance import edit_distance
from pyclustering.cluster.kmedoids import kmedoids

sequences = pd.read_csv('~/gitcode/project_amp/kmedoids/example.csv') ## read the sequence file here

## compute the pair-wise distance matrix 
num_sequences = len(sequences)
distance_matrix = np.zeros([num_sequences, num_sequences])
for i in range(num_sequences):
    for j in range (i, num_sequences):
        distance_matrix[i,j]= edit_distance(sequences['sequence'][i], sequences['sequence'][j])
distance_matrix = distance_matrix + distance_matrix.T

## compute the medoids ##
num_medoids = 10
init_medoids =  np.random.randint(num_sequences, size=num_medoids) # choose initial medoids
km = kmedoids(distance_matrix, init_medoids, data_type = 'distance_matrix')
km.process()
centers = km.get_medoids() # gets the id of the samples
clusters = km.get_clusters()

print(sequences['sequence'][centers]) # fetch the sequences which are medoids

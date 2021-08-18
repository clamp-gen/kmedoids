# https://pyclustering.github.io/docs/0.10.1/html/d0/dd3/classpyclustering_1_1cluster_1_1kmedoids_1_1kmedoids.html

import time
import os
from nltk import edit_distance
from pyclustering.cluster.kmedoids import kmedoids
import numpy as np
import pandas as pd

np.random.seed(1234)


def distance_measuring_(seq1, seq2):
    edit_value = edit_distance(seq1, seq2)
    return edit_value


def get_key(x):
    return x[1]


def _get_init_centers(n_clusters, samples):
    """return random points as initial centers"""
    # randomly choose one initial sequence as the first center
    init_ids = [np.random.randint(0, len(samples))]
    init_sequence = samples[init_ids[0]]  # sequence
    all_distances = []
    # the idea is to choose medoids far away from each other
    for _ in range(1, len(samples)):
        current_seq = samples[_]
        all_distances.append((_, distance_measuring_(init_sequence, current_seq)))
    arranged = sorted(all_distances, key=get_key,
                      reverse=True)  # descending order
    arranged_ids = [_[0] for _ in arranged]
    init_ids += arranged_ids[:n_clusters - 1]
    return init_ids


def get_distance_matrix(samples):
    dist_mat = np.zeros((len(samples), len(samples)))
    start = time.time()
    for j in range(len(samples)):
        for i in range(len(samples)):
            if i == j:
                dist_mat[i, j] = 0.
            else:
                dist_mat[i, j] = distance_measuring_(samples[i], samples[j])
    end = time.time()
    delay = (end - start) / 60
    print('Time of execution: {} minutes'.format(delay))
    return dist_mat


dataset = pd.read_csv('D1.csv')
sequences = dataset.Sequences.tolist()

initial_medoids = _get_init_centers(500, sequences)
if not os.path.exists('distance_matrix.npy'):
    print('Building the matrix distance...')
    distance_matrix = get_distance_matrix(sequences)
    np.save('distance_matrix', distance_matrix)
else:
    print('Existing distance matrix --- loading...')

def_distance_matrix = np.load('distance_matrix.npy')
matrix = def_distance_matrix.tolist()

start = time.time()
km = kmedoids(matrix, initial_medoids, data_type='distance_matrix')
km.process()
centers = km.get_medoids()
clusters = km.get_clusters()
end = time.time()
delay = (end - start) / 60
print('Time of K-medoid (PAM) execution: {} minutes'.format(delay))

print(centers)
print(clusters)

np.save('medoids', centers)  # saves final medoids ids --- we can easily retrieve the sequences later using those ids
np.save('clusters', clusters)  # saves final clusters members ids

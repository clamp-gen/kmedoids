# https://pyclustering.github.io/docs/0.10.1/html/d0/dd3/classpyclustering_1_1cluster_1_1kmedoids_1_1kmedoids.html

import time
import os
from nltk import edit_distance
from pyclustering.cluster.kmedoids import kmedoids
import numpy as np
import pandas as pd
import argparse

np.random.seed(1234)

ap = argparse.ArgumentParser()
ap.add_argument("-k", "--n_medoids", default=500, type=int,
                help="number of medoids (prototypes)")
args = vars(ap.parse_args())


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


def main(args):
    dataset = pd.read_csv('D1.csv')
    amps = dataset[dataset['Label'] == 1]['Sequences']
    sequences = amps.tolist()
    k = args.n_medoids
    initial_medoids = _get_init_centers(k, sequences)
    if not os.path.exists('distance_matrix_D1.npy'):
        print('Building the matrix distance...')
        distance_matrix = get_distance_matrix(sequences)
        np.save('distance_matrix_D1', distance_matrix)
    else:
        print('Loading existing distance matrix...')

    if not os.path.exists('medoids_k_{}.npy'.format(k)):
        def_distance_matrix = np.load('distance_matrix_D1.npy')
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

        np.save('medoids_k_{}'.format(k),
                centers)  # saves final medoids ids --- we can easily retrieve the sequences later using those ids
        np.save('clusters_k_{}'.format(k), clusters)  # saves final clusters members ids

        sequences_medoids = dataset["Sequences"][centers].tolist()
        print(sequences_medoids)

        # filtering out medoids who appears to be members of their own clusters
        clusters_sequences = {med: dataset["Sequences"][list(set(members).difference(set(centers)))].tolist() for
                              med, members in zip(sequences_medoids,
                                                  clusters)}
        print(clusters_sequences)

    else:
        print('Loading existing k-medoids ...')
        medoids = np.load('medoids_k_{}.npy'.format(k), allow_pickle=True)
        clusters = np.load('clusters_k_{}.npy'.format(k), allow_pickle=True)

        sequences_medoids = dataset["Sequences"][medoids].tolist()
        print(sequences_medoids)

        # filtering out medoids who appears to be members of their own clusters
        clusters_sequences = {med: dataset["Sequences"][list(set(members).difference(set(medoids)))].tolist() for
                              med, members in zip(sequences_medoids,
                                                  clusters.tolist())}
        print(clusters_sequences)


if __name__ == '__main__':
    args = ap.parse_args()
    main(args)

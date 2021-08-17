import numpy as np
from copy import deepcopy
import pandas as pd
from nltk.metrics.distance import edit_distance

# Just a naive way to design unique (x, y) for the peptides sequences
sequences = pd.read_csv('D1.csv')
sequences.drop_duplicates(keep=False, inplace=True)
sequences = sequences.Sequences.tolist()

uniq_char = [c for sequence in sequences for c in sequence]
char_set = {}
int_map = 1
for c in uniq_char:
    if c not in list(char_set.keys()):
        char_set[c] = int_map
        int_map += 1


def convert_sequence_to_tuple(seq_, mapper):
    # x will be overall average of the product of the mapper * the index of each character
    # in that way we create unique integer for each sequence
    # y will be overall average of the product of the mapper * the count of each character
    # a pair of (x, y) is hence unique
    x = np.array([mapper[char] * (seq_.index(char) + 1) for char in seq_])  # +1 is added to avoid index 0 effect
    x = x.mean()
    y = np.array([mapper[char] * (seq_.count(char)) for char in seq_])
    y = y.mean()
    return x, y


frame = pd.DataFrame()
xs, ys = [], []
for seq in range(len(sequences)):
    x_, y_ = convert_sequence_to_tuple(sequences[seq], char_set)
    xs.append(x_)
    ys.append(y_)

frame['seqs'] = sequences
frame['x'] = xs
frame['y'] = ys
frame.drop_duplicates(subset=["x", "y"], inplace=True)  # drops 68 sequences out of 4262, remaining has unique (x,y)
frame.to_csv('seq_x_y_for_clustering.csv')

# =============================================================== K-medoids ===============================================================================


def distance_measuring(seq1, seq2):
    edit_value = edit_distance(seq1, seq2)
    return edit_value


def get_key(x):
    return x[1]


def _get_init_centers(n_clusters, samples):
    """return random points as initial centers"""
    # randomly choose one initial sequence as the first center
    init_ids = [0]
    init_sequence = samples[0][0]  # sequence
    all_distances = []
    # the idea is to choose medoids far away from each other
    for _ in range(1, len(samples)):
        current_seq = samples[_][0]
        all_distances.append((_, distance_measuring(init_sequence, current_seq)))
    arranged = sorted(all_distances, key=get_key, reverse=True)  # descending order
    arranged_ids = [_[0] for _ in arranged]
    init_ids += arranged_ids[:n_clusters - 1]
    return init_ids


def _get_distance(data1, data2):
    """edit distance"""
    data1_seq, data2_seq = data1[0], data2[0]
    return distance_measuring(data1_seq, data2_seq)


def _get_cost(X, centers_id, dist_func):
    """return total cost and cost of each cluster"""
    dist_mat = np.zeros((len(X), len(centers_id)))
    # compute distance matrix
    for j in range(len(centers_id)):
        center = X[centers_id[j]][0]  # sequence X[centers_id[j], :]
        for i in range(len(X)):
            if i == centers_id[j]:
                dist_mat[i, j] = 0.
            else:
                dist_mat[i, j] = dist_func(X[i][0], center)  # dist_func(X[i, :], center)
    mask = np.argmin(dist_mat, axis=1)
    members = np.zeros(len(X))
    costs = np.zeros(len(centers_id))
    for i in range(len(centers_id)):
        mem_id = np.where(mask == i)
        members[mem_id] = i
        costs[i] = np.sum(dist_mat[mem_id, i])
    return members, costs, np.sum(costs), dist_mat


def _kmedoids_run(X, n_clusters, dist_func, max_iter=1000, verbose=True):
    """run algorithm return centers, members, and etc."""
    # Get initial centers
    n_samples, n_features = len(X), 2
    init_ids = _get_init_centers(n_clusters, X)
    if verbose:
        print('Initial centers are: {}'.format(init_ids))
    centers = init_ids
    members, costs, tot_cost, dist_mat = _get_cost(X, init_ids, dist_func)
    cc, swapped = 0, True
    while True:
        swapped = False
        for i in range(n_samples):
            if i not in centers:
                for j in range(len(centers)):
                    centers_ = deepcopy(centers)
                    centers_[j] = i
                    members_, costs_, tot_cost_, dist_mat_ = _get_cost(X, centers_, dist_func)
                    if tot_cost_ < tot_cost:
                        members, costs, tot_cost, dist_mat = members_, costs_, tot_cost_, dist_mat_
                        centers = centers_
                        swapped = True
                        if verbose:
                            print('Change centers to {}'.format(centers))
        if cc > max_iter:
            if verbose:
                print('End Searching by reaching maximum iteration')
            break
        if not swapped:
            if verbose:
                print('End Searching by no swaps')
            break
        cc += 1
    return centers, members, costs, tot_cost, dist_mat


class KMedoids(object):
    def __init__(self, dataset, n_clusters, dist_func=_get_distance, max_iter=300):
        self.dataset = dataset
        self.n_clusters = n_clusters
        self.dist_func = dist_func
        self.max_iter = max_iter
        self.outputs = {}
        self.medoids = {}
        self.centers_ids = []
        self.medoids_sequences = {}

    def get_key_by_value(self, value):
        seq = ''
        for _, x, y in self.dataset:
            if [x, y] == value:
                seq = _
                break
        return seq

    def fit(self, verbose=True):
        centers, members, costs, tot_cost, dist_mat = _kmedoids_run(
            self.dataset, self.n_clusters, self.dist_func, max_iter=self.max_iter, verbose=verbose)
        self.centers_ids = centers

        x = []
        for i in self.dataset:
            x.append([i[1], i[2]])
        x = np.array(x)
        for i in range(len(centers)):
            X_c = x[members == i, :]
            self.outputs[centers[i]] = X_c
            cluster_sequences = []
            temp = [X[centers[i]][1], X[centers[i]][2]]
            for seq_x_y in X_c.tolist():
                if seq_x_y != temp:
                    cluster_sequences.append(self.get_key_by_value(seq_x_y))
            self.medoids[centers[i]] = temp
            self.medoids_sequences[X[centers[i]][0]] = list(set(cluster_sequences))

        return self.outputs, self.medoids, self.medoids_sequences

    # to compute the optimal number of clusters to use
    def wcss(self):
        wcss = 0
        for k in range(self.n_clusters):
            medoids_x_y = np.array(self.medoids.values())
            wcss += np.sum((self.outputs[self.centers_ids[k]] - medoids_x_y) ** 2)
        return wcss


set_ = pd.read_csv('seq_x_y_for_clustering.csv')
X = [(seq, x, y) for seq, x, y in zip(set_.seqs.tolist(), set_.x.tolist(), set_.y.tolist())]

K = 800
model = KMedoids(X, n_clusters=K)
medoids_members, medoids_centers, seqs = model.fit(verbose=True)
np.save('clusters_medoids_sequences_and_neighbours', np.array(seqs))

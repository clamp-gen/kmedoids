import numpy as np
from nltk import edit_distance


def mode_entropy(d_hard, batch, metric=edit_distance):
    """
    :param d_hard: list of prototypes
    :param metric: distance metric --- default is edit distance
    :param batch: batch of sequences (list)
    :return: return the entropy (diversity of coverage) of each hard sequence
    """
    counts = np.zeros(len(d_hard))
    for s in batch:
        distances = [metric(s_i, s) for s_i in d_hard]
        closest = np.argmin(distances)
        counts[closest] += 1
    n = len(batch)
    p_hat = counts / n
    entropy = (p_hat * np.log(p_hat)).sum()
    return entropy


def mode_distance_average(batch, metric=edit_distance):
    """
    :param batch: batch of sequences (list)
    :param metric: distance metric - default is edit distance
    :return: average of all pairwise distances across the batch
    """
    distances = []
    for i in range(len(batch)):
        for j in range(i + 1, len(batch)):
            distances.append(metric(batch[i], batch[j]))
    return np.array(distances).mean()


def mode_count():
    """
    To be implemented
    :return:
    """

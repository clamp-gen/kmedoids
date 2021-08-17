from nltk.metrics.distance import edit_distance
import numpy as np
import pandas as pd


def o_hard_sequence(s, d_hard, c=35):
    """
    Return O_hard similarity to d_hard set
    :param s: sequence to evaluate
    :param d_hard: set of k medoids
    :param c: constant of normalisation
    :return: the score characterizing how close that sequence is to all main k medoids
    """
    distances = [edit_distance(s, s_) for s_ in d_hard]
    exp_arg = min(distances) / c
    return np.exp(-exp_arg)


def o_hard_batch(batch, d_hard, top_k=1000):
    """
    Given a batch of n sequences, return the top k sequences based on o_hard scores
    :param batch: list of sequences
    :param d_hard: set of k medoids
    :param top_k: number of sequences to retain after the filtering
    :return: dataframe of top_k sequences and their respective scores
    """
    filtered = pd.DataFrame()
    seqs, scores = [], []
    for s in batch:
        seqs.append(s)
        scores.append(o_hard_sequence(s, d_hard))
    filtered['seqs'] = seqs
    filtered['o_hard_scores'] = scores
    filtered = filtered.sort_values(by=['o_hard_scores'], ascending=False)
    return filtered[:top_k]

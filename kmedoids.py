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

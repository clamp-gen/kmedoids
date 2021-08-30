from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from diamond_command import DiamondCommandline
import time
import numpy as np
import pandas as pd


def get_distance_matrix(samples):
    dist_mat = np.zeros((len(samples), len(samples)))
    start = time.time()
    for i in range(len(samples)):
        for j in range(i, len(samples)):
            if i == j:
                dist_mat[i, j] = 0.
            else:
                dist_mat[i, j] = diamond_e_value(samples[i], samples[j])
    end = time.time()
    delay = (end - start) / 60
    print('Time of execution: {} minutes'.format(delay))
    return dist_mat


def diamond_e_value(seq1, seq2):
    """ returns the e-value of two sequences using BLAST"""
    # The Expect value (E) is a parameter that describes the number of hits one can "expect" to see by chance when
    # searching a database of a particular size. It decreases exponentially as the Score (S) of the match increases.
    # Thus E value is a measure of distance between two sequences

    # Create two sequence files
    seq1 = SeqRecord(Seq(seq1),
                     id="seq1")
    seq2 = SeqRecord(Seq(seq2),
                     id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")
    db_name = "seq2"

    DiamondCommandline(createdb=True, infile="seq2.fasta", db=db_name)()

    # Run BLAST and parse the output as XML

    output = DiamondCommandline(bp=True, query="seq1.fasta", db=db_name, outfmt=6)()[0]

    # TODO: develop a XML based read
    # outfmt 6 is a tabular with default format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    if output:
        evalue = float(str(output).split('\t')[10])
    else:
        evalue = 1e100  # set a very large value

    return evalue

# Test case #
def main():
    dataset = pd.read_csv('D1.csv')
    sequences = dataset.Sequences.tolist()
    amp = dataset.loc[dataset['Label'] == 1]
    amp = amp.Sequences.tolist()
    # amp = amp[0:100]
    start = time.time()
    dist_mat = get_distance_matrix(amp)
    end = time.time()
    print(" time use ", end - start)
    np.save('diamond_distance_matrix', dist_mat)


if __name__ == "__main__":
    main()
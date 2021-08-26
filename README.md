# kmedoids

Precomputed distance metric on our current set D1 (used for the real oracle) accessible [here](https://drive.google.com/file/d/17I8He-Dn3vXu3q3b7Ge_h5bw9GJ5yTgC/view?usp=sharing)


To run: python k_medoids_better.py --n_medoids 500

- This will load the default distance matrix built on the full D1 set (real oracle) - the list of medoids (initially 500), and the full list of associated neighbors (clusters members)
  - If n_medoids != 500, the distance matrix will be used and new set of k-medoids will be created
  - The 500 k medoids (protoypes) are available [here](https://github.com/clamp-gen/kmedoids/blob/main/medoids_k_500.npy)
  - The set of members of the cluster of each medoid is available [here](https://github.com/clamp-gen/kmedoids/blob/main/clusters_k_500.npy)

The HaackMD note summarizing the approach is [here](https://hackmd.io/hoWSQvl5Tl6CsxoPuANKyw)

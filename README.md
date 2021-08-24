# kmedoids

Precomputed distance metric on our current set D1 (used for the real oracle) accessible [here](https://drive.google.com/file/d/12zqAvnPKjUBvbP0O405CkL0-zY-iLudq/view?usp=sharing)


To run: python k_medoids_better.py --n_medoids 500

- This will load the default distance matrix built on the full D1 set (real oracle) - the list of medoids (initially 500), and the full list of associated neighbors (clusters members)
  - If n_medoids != 500, the distance matrix will be used and new set of k-medoids will be created

The HaackMD note summarizing the approach is [here](https://hackmd.io/hoWSQvl5Tl6CsxoPuANKyw)

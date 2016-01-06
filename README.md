### About

This is a really barbones toy project to illustrate how to compute a Gram matrix for the one-sided mean alignment kernel.

It reads a set of sequences from file `example_data.txt`, computes a Gram matrix of these sequences using the _one-sided mean alignment kernel_, computes the eigenvalues of this matrix and displays them.
As the one-sided mean alignment kernel is provably positive definite, the eigenvalues are all positive.

The project also contains code to generate random dummy time series.

### Usage

The project doesn't need any dependency so you can just execute it:

```bash
python mean_alignment.py
```

### About the one-sided mean alignment kernel

Well-known kernel methods such as for example SVM and Kernel PCA rely on a kernel which is used computes a Gram matrix which can generally be interpreted as a matrix of similarity values between any pair of samples in the dataset.
These algorithms require that the kernel be _positive definite_, which guarantees that the resulting optimization programs are convex whatever the samples.

When dealing with vector data the most sensible choice is generally a Gaussian kernel, but when the samples are time-series (or more generally sequences) custom kernels must be used.
There are not many sensible choices, as merely using a classic _dynamic time warping_ distance does not lead to a positive definite kernel.

To this end we propose the one-sided mean kernel, which has many advantages:
* Provably positive definite,
* Faster than competing techniques with a time complexity of `O(l × (m - l))` instead of `O(l × m)` for sequences of length `l < m`,
* Consistent with a vector kernel when time series are of equal length,
* Does not suffer from issues of diagonal dominance like the global alignment kernel for example.

An implementation in Haskell is available [here](https://github.com/nchrys/haskell-sequence-kernels).

### References

* N. Chrysanthos, P. Beauseroy, H. Snoussi, E. Grall-Maës, __Theoretical properties and implementation of the one-sided mean kernel for time series__, in: _Neurocomputing_ 169 (2015) 196–204
* M. Cuturi, J.-P. Vert, O. Birkenes, T. Matsui, __A kernel for time series based on global alignments__, in: _IEEE International Conference on Acoustics, Speech and Signal Processing_, 2007
* J.-P. Vert, __The Optimal Assignment Kernel is not Positive Definite__ in: _arXiv:0801. 4061_

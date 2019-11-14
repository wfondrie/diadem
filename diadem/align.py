"""
This module contains the my implementation of the FastDTW algorithm.

The algorithm is described in http://cs.fit.edu/~pkc/papers/tdm04.pdf.
This implementation is losely based on the python package from this
GitHub repository: https://github.com/slaypni/fastdtw.

My code deviates from this repository is a few ways to make it more
user friendly and amenable to aligning mass spectrometry runs:
  1. Cython is not needed for good speed, because of numba.
  2. The input numpy arrays (x and y) can be of any dimensionality, so
     long as the distance function can handle it.

Written by William E Fondrie, 2019
"""
from typing import Tuple, Callable

import numpy as np
import numba as nb

# Distance Functions ----------------------------------------------------------
@nb.njit
def cosine_distance(x, y, tiny=np.finfo(float).tiny):
    """Compute 1 minus the cosine similarity between x and y"""
    denom = (np.linalg.norm(x) * np.linalg.norm(y) + tiny)
    return 1 - np.dot(x, y) / denom


# DTW Functions ---------------------------------------------------------------
def fastdtw(x: np.ndarray, y: np.ndarray, radius: int = 1,
            dist: Callable[[np.ndarray, np.ndarray], float]
            = cosine_distance) -> Tuple[float, Tuple[Tuple[int, int]]]:
    """
    Find the approximate minimum warping path between x and y.

    Parameters
    ----------
    x, y : numpy.ndarray
        Numpy arrays of the series to align. The first dimension is
        always assumed to be the time domain. For example, if aligning
        two mass spectrometry runs by their precursor mass spectra,
        x and y would be of shape [retention time, m/z] where m/z is
        each spectrum vectorized along the m/z axis.

    radius : int
        The radius to use for the FastDTW neighborhood.

    dist: Callable
        A distance function (not in the strict sense of the word), which
        accepts single time slices of x and y as input and returns their
        distance as a float.

    Returns
    -------
    Tuple[float, Tuple[Tuple[int, int]]]
        A tuple containing two elements. The first is the estimated DTW
        distance between x and y. The second is a Tuple of Tuples
        indicating the minimal warping path between x and y. The
        innermost tuple contains the mapping of (x, y) pairs in the
        path.
    """
    min_time_size = radius + 2

    # The base case
    if x.shape[0] < min_time_size or y.shape[0] < min_time_size:
        return dtw(x, y, dist)

    # Recursive state
    shrunk_x = _reduce_by_half(x)
    shrunk_y = _reduce_by_half(y)
    _, path = fastdtw(shrunk_x, shrunk_y, radius=radius)
    window = _expand_window(path, x.shape[0], y.shape[0], radius=radius)

    return dtw(x, y, dist, window)


def dtw(x: np.ndarray, y: np.ndarray,
        dist: Callable[[np.ndarray, np.ndarray], float] = cosine_distance,
        _window = None) -> Tuple[float, Tuple[Tuple[int, int]]]:
    """
    Find the minimum warping path between x and y.

    Parameters
    ----------
    x, y : numpy.ndarray
        Numpy arrays of the series to align. The first dimension is
        always assumed to be the time domain. For example, if aligning
        two mass spectrometry runs by their precursor mass spectra,
        x and y would be of shape [retention time, m/z] where m/z is
        each spectrum vectorized along the m/z axis.

    dist: Callable
        A distance function (not in the strict sense of the word), which
        accepts single time slices of x and y as input and returns their
        distance as a float.

    Returns
    -------
    Tuple[float, Tuple[Tuple[int, int]]]
        A tuple containing two elements. The first is the estimated DTW
        distance between x and y. The second is a Tuple of Tuples
        indicating the minimal warping path between x and y. The
        innermost tuple contains the mapping of (x, y) pairs in the
        path.
    """
    if _window is None:
        _window = [(i, j) for i in range(x.shape[0]) for j in range(y.shape[0])]

    _window = list(_window)
    return _dtw_main(x, y, dist, _window)

# Utility functions -----------------------------------------------------------
# This is the implementation of the Dynamic Time Warping algorithm.
# For some reason the jitted version is wayyyyy slower :(
def _dtw_main(x, y, dist, window):
    """The DTW algorithm"""
    res = {}
    res[0, 0] = (float(0), 0, 0)

    for i, j in window:
        dt = dist(x[i, ...], y[j, ...])
        moves = ((i, j+1), (i+1, j), (i, j))

        val = np.Inf
        for move in moves:
            if move in res:
                if res[move][0] < val:
                    val = res[move][0]
                    res[i+1, j+1] = (val + dt, *move)


    path = []
    i, j = x.shape[0], y.shape[0]
    while i or j:
        path.append((i-1, j-1))
        i, j = res[i, j][1], res[i, j][2]

    path.reverse()
    return (res[x.shape[0], y.shape[0]][0], tuple(path))


def _reduce_by_half(x):
    """Reduce x by half by taking the average."""
    max_idx = x.shape[0] - (x.shape[0] % 2)
    return np.array([(x[i, ...] + x[i+1, ...]) / 2 for i in range(0, max_idx, 2)])


def _expand_window(path, len_x, len_y, radius):
    """Expands the window around path and returns a new window"""
    path_ = set(path)
    path_range = range(-radius, radius+1)
    window = set()

    for i, j in path:
        for a, b in ((i+a, j+b) for a in path_range for b in path_range):
            if 0 <= a < len_x and 0 <= b < len_y:
                path_.add((a, b))

    for i, j in path_:
        i *= 2
        j *= 2
        for a, b in ((i, j), (i, j+1), (i+1, j), (i+1, j+1)):
            if 0 <= a < len_x and 0 <= b < len_y:
                window.add((a, b))

    return sorted(window)

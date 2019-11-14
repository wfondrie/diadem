"""
Read mzML files.
"""
from typing import Tuple

import tqdm
import numpy as np
from pyteomics import mzml

from diadem.dataset import DIAScan, DIARun

def read(mzml_file: str, max_peaks: int = None, min_intensity: float = None)\
         -> Tuple[np.ndarray, DIAScan]:
    """
    Read an mzML file from a DIA experiment.

    Parameters
    ----------
    mzml_file : str
        The mzML file to read.

    Returns
    -------
    diadem.dataset.DIARun
        A DIARun object containg the raw data.
    """
    kwargs = {"max_peaks": max_peaks, "min_intensity": min_intensity}
    with mzml.MzML(mzml_file) as mz_dat:
        scans = DIARun([s for s in _pbar(mz_dat.map(_mkscan, kwargs=kwargs, processes=4))])

    return scans


def _pbar(x):
    """Create a tqdm progress bar"""
    return tqdm.tqdm(x, ascii=True, unit=" scans")


def _mkscan(spectrum, max_peaks, min_intensity):
    """Read a scan"""
    scan = DIAScan(spectrum)
    scan.preprocess(min_intensity, max_peaks)
    return scan

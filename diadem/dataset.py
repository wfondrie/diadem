"""
Classes for storing and manipulating DIA datasets.

The DIARun class stores data from a single mzML file from a DIA
experiment. It consists of a collection of DIAWindow() instances which
store data from the individual DIA MS/MS windows.
"""
import typing

import numpy as np
import pandas as pd

class DIARun():
    """
    Store and manipulate a DIA run.

    Initialization loads a data-independent acquisition mass
    spectrometry (DIA) mzML file into memory.

    Parameters
    ----------
    mzml_file : str
        The mzML file to read. This can be either a normal mzML file or
        gzipped (mzML.gz).


    Attributes
    ----------
    windows : 

    """
    def __init__(self):
        pass


class DIAWindow():
    """
    Store data for each DIA MS/MS isolation window.

    Parameters
    ----------
    mzml_data : dict
        A dictionary with the required data from the mzML file.

    Attributes
    ----------
    name : str
        Name for the DIAWindow.

    index : numpy.ndarray
        A 1D numpy array of scan numbers.

    scan : Tuple[str]
        The unique identifier for each scan.

    peaks : Tuple[numpy.ndarray]
        A tuple of 1D numpy arrays, containing the mass peaks of each
        scan.

    ret_time : numpy.ndarray
        A 1D numpy array with the raw retention time of each scan.

    tic : np.ndarray
        A 1D numpy array with the tic of each scan.
    """
    def __init__(self, name, mz_mid, mz_width):
        pass


class DIAScan():
    """
    Store data for an individual DIA Scan.
    """
    def __init__(self, spectrum):
        """Initialize a DIAScan() object"""
        ms_level = spectrum["ms level"]

        if ms_level == 1:
            self.window = "precursor"

        else:
            mz_info = spectrum["precursorList"]["precursor"]
            mz_info = mz_info[0]["isolationWindow"]
            mz_info = mz_info["isolation window target m/z"]
            self.window = f"{mid:.0f}"


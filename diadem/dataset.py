"""
Classes for storing and manipulating DIA datasets.

The DIARun class stores data from a single mzML file from a DIA
experiment. It consists of a collection of DIAWindow() instances which
store data from the individual DIA MS/MS windows.
"""
import logging
import multiprocessing as mp
from typing import Tuple, Dict
from itertools import chain

import tqdm
import numpy as np
import numba as nb

import diadem.write
from diadem.align import fastdtw

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
    def __init__(self, scans, spectrum_matrix=None):
        """Initialize a DIARun"""
        self.spectrum_matrix = spectrum_matrix
        self.windows = {}
        for scan in scans:
            win = scan.window[0]
            curr_win = self.windows.get(win, [])
            curr_win.append(scan)
            self.windows[win] = curr_win

        for window, spec_list in self.windows.items():
            self.windows[window] = DIAWindow(window, spec_list)

    @property
    def scans(self):
        """Get the scans, ordered by index"""
        all_scans = [s.scans for _, s in self.windows.items()]
        all_scans = list(chain.from_iterable(all_scans))
        all_scans.sort(key=_get_index)
        return all_scans

    def mask(self, reference_run, tol=10):
        """Mask this run using a reference run"""
        logging.info("Masking run by window...")
        num_filtered = []
        for win_name, win in tqdm.tqdm(self.windows.items()):
            ref_win = reference_run.windows[win_name]
            before, after = win.mask(ref_win, tol)
            num_filtered.append((win_name, before, after))

        total_before = 0
        total_after = 0
        for win in num_filtered:
            logging.info("%s: %i -> %i peaks (%0.2f%% remaining)",
                         win[0], win[1], win[2], win[2]/win[1]*100)
            total_before += win[1]
            total_after += win[2]

        logging.info("Total: %i -> %i peaks (%0.2f%% remaining)",
                     total_before, total_after, total_after/total_before*100)

    def align(self, reference_run, radius=1):
        """Calibrate the retention time to a reference run."""
        logging.info("Aligning runs by window...")
        for win_name, win in tqdm.tqdm(self.windows.items()):
            ref_mat = reference_run.windows[win_name].vectorize()
            targ_mat = win.vectorize()
            _, path = fastdtw(ref_mat, targ_mat, radius=radius)
            indices = _path2map(path, len(win.scans))
            win.reference_index = indices

    def write(self, out_file):
        """Write a DIARun to mzML"""
        diadem.write.mzML(self, out_file)


class DIAWindow():
    """
    Store data for each DIA MS/MS isolation window.

    Parameters
    ----------
    scans : A list of DIAScan
    """
    def __init__(self, name, scans):
        self.name = name
        self._scans = scans
        self.sort_scans()

    @property
    def tic(self):
        """The tic of each scan"""
        return np.array([s.tic for s in self._scans])

    @property
    def original_tic(self):
        return np.array([s.original_tic for s in self._scans])

    @property
    def peaks(self):
        return np.array([s.peaks for s in self._scans])

    @property
    def original_peaks(self):
        return np.array([s.original_peaks for s in self._scans])

    @property
    def ret_time(self):
        """The retention times of each scan"""
        return np.array([s.ret_time for s in self._scans])

    @property
    def reference_index(self):
        """The reference scan index for each scan"""
        return np.array([s.reference_index for s in self._scans])

    @reference_index.setter
    def reference_index(self, indices):
        """Set the calibrated retention time for all scans"""
        for idx, scan in zip(indices, self._scans):
            scan.reference_index = idx

    @property
    def scans(self):
        """Return the scans in a window"""
        return self._scans

    def sort_scans(self):
        """Sort the scans by raw_rt"""
        self._scans.sort(key=_get_index)

    def vectorize(self, **kwargs) -> np.ndarray:
        """Vectorize the spectra and return a matrix"""
        vecs = [s.vectorize(**kwargs) for s in self.scans]
        return np.vstack(vecs)

    def mask(self, reference_window, tol):
        """Mask using a reference window"""
        total = 0
        kept = 0
        for scan in self.scans:
            ref_mz = [reference_window.scans[i].mz for i in scan.reference_index]
            before, after = scan.mask(np.concatenate(ref_mz), tol)
            total += before
            kept += after

        return (total, kept)


class DIAScan():
    """
    Store data for an individual DIA Scan.
    """
    def __init__(self, spectrum):
        """Initialize a DIAScan() object"""
        keys = ["ms level", "scanList", "m/z array",
                "intensity array", "id", "index"]

        if spectrum["ms level"] == 2:
            keys += ["precursorList"]

        self._spectrum = {key: spectrum[key] for key in keys}
        self.reference_index = [None]
        self.original_peaks = len(self._spectrum["intensity array"])
        self.original_tic = self._spectrum["intensity array"].sum()

    @property
    def ms_level(self) -> int:
        """Retrieve the MS level"""
        return self._spectrum["ms level"]

    @property
    def ret_time(self) -> float:
        """Get the retention time of the spectrum"""
        return self._spectrum["scanList"]["scan"][0]["scan start time"]

    @property
    def window(self) -> Tuple[str, float, float, float]:
        """Retrieve the DIA window."""
        if self.ms_level == 1:
            return ("precursor",)

        mz_info = self._spectrum["precursorList"]["precursor"][0]
        mz_info = mz_info["isolationWindow"]
        mid = mz_info["isolation window target m/z"]
        offset_low = mz_info["isolation window lower offset"]
        offset_high = mz_info["isolation window upper offset"]

        return (f"m/z {mid-offset_low:0.0f}-{mid+offset_high:0.0f}",
                mid, offset_low, offset_high)

    @property
    def mz(self) -> np.ndarray:
        """Retrieve the m/z values of the spectrum"""
        return self._spectrum["m/z array"]

    @mz.setter
    def mz(self, mz_array):
        """Set the m/z array"""
        self._spectrum["m/z array"] = mz_array

    @property
    def intensity(self) -> np.ndarray:
        """Retrieve the intensity values of the spectrum"""
        return self._spectrum["intensity array"]

    @intensity.setter
    def intensity(self, intensity_array):
        """Set the intensity array"""
        self._spectrum["intensity array"] = intensity_array

    @property
    def peaks(self) -> int:
        """Get the current number of peaks in the spectrum"""
        return len(self.intensity)

    @property
    def activation(self) -> Dict:
        """Retrieve the activation method for the spectrum"""
        if self.ms_level == 1:
            return None

        mz_info = self._spectrum["precursorList"]["precursor"][0]
        return mz_info["activation"]

    @property
    def scan(self) -> str:
        """Retrieve the scan header"""
        return self._spectrum["id"]

    @property
    def index(self) -> int:
        """Retrieve the scan index number"""
        return self._spectrum["index"]

    @property
    def tic(self) -> float:
        """Retrieve the total ion current of the spectrum"""
        return self.intensity.sum()

    def vectorize(self, bin_width: float = 1.0005,
                  min_mz: float = 0.4, max_mz: float = 2000.0) \
        -> np.ndarray:
        """
        Vectorize the mass spectrum

        Parameters
        ----------
        bin_width : float
            The bin width to use for vectorization
        min_mz : float
            The lowest m/z bin.
        max_mz : float
            The highest m/z bin.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of the vectorize spectrum.
        """
        return _vectorize(self.mz, self.intensity, bin_width,
                          min_mz, max_mz)
        

    def mask(self, mask_mz, tol=10) -> Tuple[np.ndarray]:
        """
        Mask the mass spectrum by removing peaks within the tolerance any peaks in the list.
        """
        # Returns the indexes of ions to keep.
        idx = _mask(self.mz, mask_mz, tol)
        before = len(self.mz)
        self.filter(idx)
        after = len(self.mz)

        return (before, after)


    def filter(self, index: np.ndarray) -> None:
        """
        Filter the m/z and intensity arrays jointly, keeping those
        specified by index

        Parameters
        ----------
        index : numpy.ndarray
            The indices of elements to keep.
        """
        self.mz = self.mz[index]
        self.intensity = self.intensity[index]


    def preprocess(self, min_intensity: float = None,
                   max_peaks: int = None) -> None:
        """
        Preprocess a mass spectrum.

        Note that this will modify the m/z and intensity arrays.

        Parameters
        ----------
        min_intensity : float
            Specify the minimal fraction the base peak intensity that a peak
            must have to be kept. None keeps all.
        max_peaks : float
            The maximum number of most intense peaks to keep. None keeps all.

        Returns
        -------
        DIAScan
            A DIAScan object that has been preprocessed.
        """
        if min_intensity is not None:
            frac = self.intensity / self.intensity.max()
            self.filter(np.nonzero(frac >= min_intensity))

        if max_peaks is not None and max_peaks < len(self.intensity):
            n = -1*max_peaks
            self.filter(np.argpartition(self.intensity, n)[n:])

# Utility Functions -----------------------------------------------------------
def _align(ref_mat, targ_mat, radius):
    """Align a window"""
    _, path = fastdtw(ref_mat, targ_mat, radius)
    return _path2map(path, len(target_window.scans))

def _get_index(scan):
    """Return the index of a scan"""
    return scan.index

def _path2map(path, targ_length):
    """
    Turn a path into a list of lists mapping target scans to one or more
    reference scans.
    """
    idx_map = [[] for x in range(targ_length)]
    for step in path:
        idx_map[step[1]].append(step[0])

    return idx_map

@nb.njit
def _mask(targ_mz, mask_mz, tol):
    """Filter the targ_mz array for peaks that do not have a match in mask_mz"""
    tol = tol * 1e-6
    ret_indices = []
    for idx, targ in enumerate(targ_mz):
        for mask in mask_mz:
            diff = (targ - mask) / targ
            within_tol = np.abs(diff) <= tol
            if within_tol:
                break

        if not within_tol:
            ret_indices.append(idx)

    return ret_indices

@nb.njit
def _vectorize(mz_array, intensity_array, bin_width,
               min_mz, max_mz):
    """Quickly vectorize a spectrum"""
    bins = np.arange(min_mz, max_mz, bin_width)
    bin_idx = np.digitize(mz_array, bins)
    unique_idx = np.unique(bin_idx)
    bin_int = [np.max(intensity_array[bin_idx == x]) for x in unique_idx]
    vec = np.zeros(len(bins)+1)
    vec[unique_idx] = np.array(bin_int)

    return vec[1:-1] # trim bins outside of (min_mz, max_mz)

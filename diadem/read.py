"""
Read mzML files.
"""
from time import time

from pyteomics import mzml
from dataset import DIAScan


def test(mzml_file: str):
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
    print(mzml_file)
    with mzml.MzML(mzml_file) as mzml_data:
        #scans = [s for s in mzml_data.map(DIAScan, 1)]
        scans = []
        for f in mzml_data.map():
            scans.append(f)


    return scans


if __name__ == "__main__":
    start = time()
    x = test("test.mzML")
    finish = time()

    print(f"{finish-start:.2f}")

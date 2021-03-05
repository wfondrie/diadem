"""
Write the masked files to mzML
"""
import tqdm
import logging
import numpy as np
from psims.mzml import MzMLWriter

from diadem import __version__

def mzML(run, out_file):
    """Write a DIARun to an mzML file"""
    scans = run.scans
    writer = MzMLWriter(out_file)
    with writer:
        writer.controlled_vocabularies()
        writer.file_description(["MS1 spectrum", "MSn spectrum",
                                 "centroid spectrum"])
        writer.software_list([{"id": "diadem", "version": __version__}])
        writer.instrument_configuration_list(None)

        methods = [writer.ProcessingMethod(order=1, software_reference="diadem")]
        processing = [writer.DataProcessing(methods, id="DP1")]
        writer.data_processing_list(processing)

        with writer.run(id="diadem"):
            with writer.spectrum_list(count=len(scans)):
                for scan in tqdm.tqdm(scans, ascii=True):
                    if scan.mz.size == 0:
                        scan.mz = np.array([0])
                        scan.intensity = np.array([0])

                    params = [{"ms level": scan.ms_level},
                              {"total ion current": scan.tic}]
                    pre = {"activation": scan.activation}
                    if scan.ms_level == 1:
                        pre["spectrum_reference"] = scan.scan
                        writer.write_spectrum(scan.mz, scan.intensity,
                                              id=scan.scan,
                                              centroided=True,
                                              scan_start_time=float(scan.ret_time),
                                              params=params)
                    else:
                        pre["mz"] = scan.window[1]
                        pre["charge"] = 1
                        pre["intensity"] = 0
                        pre["isolation_window_args"] = {"target": scan.window[1],
                                                        "lower": scan.window[2],
                                                        "upper": scan.window[3]}

                        writer.write_spectrum(scan.mz, scan.intensity,
                                              id=scan.scan,
                                              centroided=True,
                                              scan_start_time=float(scan.ret_time),
                                              precursor_information=pre,
                                              params=params)
        logging.info("Closing mML...")

    logging.info("Formatting mzML...")
    writer.format()
    logging.info("Validating mzML...")
    #writer.validate()

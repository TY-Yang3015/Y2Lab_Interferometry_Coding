from LocalCalibration import LocalCalibrator
import numpy as np
import sys

file='%s'%(sys.argv[1])
results = np.loadtxt(file, dtype='float64').T

results[0] = results[1]

calibrator = LocalCalibrator(results)
calibrator.calibrate(546e-9, 1000000, "filter")
calibrator.visualise(interactive=True)

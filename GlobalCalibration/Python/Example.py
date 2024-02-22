from GlobalCalibration import GlobalCalibrator
import numpy as np
import sys

file='%s'%(sys.argv[1])
results = np.loadtxt(file, dtype='float64').T

calibrator = GlobalCalibrator(results)
calibrator.calibrate(6.33579e-11, 10000, "filter")
calibrator.visualise(interactive=True)

from GlobalCalibration import GlobalCalibrator
import numpy as np
import sys

file='%s'%(sys.argv[1])
results = np.loadtxt('./quick_data.txt', dtype='float64').T

calibrator = GlobalCalibrator(results)
calibrator.calibrate(1.94e-11, 10000, "filter")
calibrator.visualise(interactive=True)

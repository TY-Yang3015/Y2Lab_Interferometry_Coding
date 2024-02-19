from CrossingPoint import CrossingPointsAnalyser
import numpy as np
import sys

file = '%s'%(sys.argv[1])
data = np.loadtxt(file).T

analyser = CrossingPointsAnalyser(data)
analyser.run()

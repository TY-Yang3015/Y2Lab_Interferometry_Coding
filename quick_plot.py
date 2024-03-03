import sys
import numpy as np
import matplotlib.pyplot as plt


file='%s'%(sys.argv[1])
results = np.loadtxt(file).T

try:
    import seaborn as sns
    sns.set_theme()
    sns.set_context("paper")
    sns.set(rc={"xtick.bottom" : True, "ytick.left" : True})
    palette = sns.color_palette('pastel')
except ImportError:
    pass

y1 = np.array(results[0])
y2 = np.array(results[1])

x=np.array(results[5])

if y1.all() != np.zeros(len(y1)).all():
    plt.figure("Detector 1")
    plt.plot(x,y1,'o-')
    plt.xlabel("Position in $\mu$steps")
    plt.ylabel("Signal 1")
    plt.savefig("./quick_plot_detector_1.png")

if y2.all() != np.zeros(len(y2)).all():
    plt.figure("Detector 2")
    plt.plot(x,y2,'o-')
    plt.xlabel("Position in $\mu$steps")
    plt.ylabel("Signal 2")
    plt.savefig("./quick_plot_detector_2.png")

plt.show()

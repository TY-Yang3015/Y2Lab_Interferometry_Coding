import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from statistics import mode

class CrossingPointsAnalyser:
    def __init__(self, data):
        self._set_plot_style()
        self.data = data
        self.x = self.data[5]
        if self.data[0].all() == np.zeros(len(self.data[0])).all():
            self.y = self.data[1]
        else:
            self.y = self.data[0]
        self.separations = None
        self.crossing_points = None
        
        sampling_frequency = int(mode(self.data[4]))
        if sampling_frequency == 16:
            self.sampling_frequency = 50
        elif sampling_frequency == 13:
            self.sampling_frequency = 200
        elif sampling_frequency == 11:
            self.sampling_frequency = 500
        else: raise ValueError("unable to identify sampling frequency")
        
        self.x = self.x[self.data[4] == sampling_frequency]
        self.y = self.y[self.data[4] == sampling_frequency]
        
    @classmethod
    def _set_plot_style(cls):
        try:
            import seaborn as sns
            sns.set_theme()
            sns.set_context("paper")
            sns.set(rc={"xtick.bottom" : True, "ytick.left" : True})
            palette = sns.color_palette('pastel')
        except ImportError:
            pass
        
    def _filter_y(self, y:np.ndarray, filter_order:int=2, freq:float=1) -> np.ndarray:
    
        return signal.sosfilt(signal.butter(filter_order, freq, 'hp', fs=self.sampling_frequency,
                                 output='sos'), y)
    
    def _remove_artefact(self, x:np.ndarray, y:np.ndarray, percent:float=0.05) -> np.ndarray:
        y = y[int(len(y)*(percent)):int(len(y)*(1-percent))]
        x = x[int(len(x)*(percent)):int(len(x)*(1-percent))]
        x -= min(x)
        return x, y
    
    def _calc_crossing_points_separation_adjusted(self, x:np.ndarray, y:np.ndarray) -> np.ndarray:
        crossing_points = np.array([])
        for i in range(len(y)-1):
            if (y[i] <= 0 and y[i+1] >= 0) or (y[i] >= 0 and y[i+1] <= 0):
                m = (y[i+1] - y[i]) / (x[i+1] - x[i])
                c = y[i] - m * x[i]
                if m != 0:
                    crossing_x = -c / m
                    crossing_points = np.append(crossing_points, crossing_x)
        
        self.crossing_points = crossing_points
        return np.abs(np.diff(crossing_points))
    
    def _calc_crossing_points_separation(self, x:np.ndarray, y:np.ndarray) -> np.ndarray:
        crossing_points = np.array([])
        for i in range(len(y)-1):
            if (y[i] <= 0 and y[i+1] >= 0) or (y[i] >= 0 and y[i+1] <= 0):
                b = (y[i+1] - y[i]/x[i] * x[i+1])/(1-x[i+1]/x[i])
                a = (y[i] - b)/x[i]
                extra = -b/a - x[i]
                crossing_points = np.append(crossing_points, (x[i]+extra))
    
        self.crossing_points = crossing_points
        return np.abs(np.diff(crossing_points))
        
    def run(self, reference:float=532e-9, preprocessing:bool=True):
        
        if preprocessing is True:
            self.y = self._filter_y(self.y)
            self.x, self.y = self._remove_artefact(self.x, self.y)
        elif preprocessing is False:
            self.y -= np.mean(self.y)
        else: raise ValueError("only booleans are accepted for 'preprocessing'.")
        

        self.separations = self._calc_crossing_points_separation_adjusted(self.x, self.y)
        print("The mean difference between crossing points is %.0d and the standard\
          deviation between crossing points is %.0d" % (np.mean(self.separations),
          np.std(self.separations, ddof=1)))
        
        mean_diff = np.mean(self.separations) # difference in musteps between crossing points
        metres_per_microstep = reference / (2.0*mean_diff)

        print('Distance moved per mu step is %.5e metres' % (metres_per_microstep))

        
    def visualise(self):
        fig, ax = plt.subplots(2, 1, figsize=(8, 7))
        fig.suptitle('Crossing Point Analysis', fontsize=20)

        ax[0].plot(self.x, self.y, 'x-', label='Data', color='blue')
        ax[0].plot(self.crossing_points, 0*np.array(self.crossing_points),
                   'o', label='Crossing Points', color='red')
        ax[0].set_xlabel("Position (µsteps)")
        ax[0].set_xlim(0, 0.025e7)
        ax[0].set_ylim(-500, 500)
        ax[0].set_ylabel("Signal")
        ax[0].set_title('Position against Signal (zoomed in)')
        ax[0].legend(loc='upper right')
        ax[1].set_title('Histogram of the Crossing Point Distances')
        ax[1].hist(self.separations, bins=100, color='blue')
        ax[1].set_xlabel("Distance between crossings (µsteps)")
        ax[1].set_ylabel("Number of entries")
        fig.tight_layout()

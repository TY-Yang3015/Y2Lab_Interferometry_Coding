import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate
import warnings
import plotly.express as px
import pandas as pd
from statistics import mode


class GlobalCalibrator():
    def __init__(self, data:np.ndarray):
        self._set_plot_style()
        self._set_warning()
        
        self.data = data.astype('float64')
        self.x = self.data[5].astype('float64')
        if self.data[0].all() == np.zeros(len(self.data[0])).all():
            self.y = self.data[1].astype('float64')
        else:
            self.y = self.data[0].astype('float64')
            
        self.remove_offset=None
       
        sampling_frequency = int(mode(self.data[4]))
        if sampling_frequency == 16:
            self.sampling_frequency = 50
        elif sampling_frequency == 13:
            self.sampling_frequency = 200
        elif sampling_frequency == 11:
            self.sampling_frequency = 500
            
        self.x = self.x[self.data[4] == sampling_frequency]
        self.y = self.y[self.data[4] == sampling_frequency]
            
        self.lambda_ = None
        self.yf = None
       
        self.metres_per_microstep = None
        
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
            
    @classmethod
    def _set_warning(cls):
        
        def custom_formatwarning(msg, category, *args, **kwargs):
            return f'{category.__name__}: {msg}\n'
        
        warnings.formatwarning = custom_formatwarning
        
    
    def _auto_sanitise(self, x0:np.ndarray, y0:np.ndarray) -> np.ndarray:
        
        dx = np.diff(x0)
        anomalies = np.argwhere((dx < 0) | (dx == 0))
        
        percentage = 0
        #decide if sanitise
        if len(anomalies) == 0:
            print("no sanitisation applied.")
            return x0, y0
            
        else:
            #decide strategy: one-tail or two-tail
            if np.std(anomalies) < (len(x0)/3):
                #one-tail: decide cut head or tail
                if np.mean(anomalies) > (len(x0)/2):
                    x = x0[:np.min(anomalies)]
                    y = y0[:np.min(anomalies)]
                    percentage = (float(len(x0)) - float(np.min(anomalies)))/float(len(x0))
                    print("one-tail sanitisation at tail applied.")
                    print("number of sanitised elements: ", int(len(x0) - np.min(anomalies)))
                    print("number of anomalous elements: ", len(anomalies))
                    print("sanitisation percentage: ", format(percentage*100., '.2e'), "%")
                else:
                    x = x0[np.max(anomalies)+1:]
                    y = y0[np.max(anomalies)+1:]
                    percentage = (np.max(anomalies))/len(x0)
                    print("one-tail sanitisation at head applied.")
                    print("number of sanitised elements: ", int(np.max(anomalies) +1))
                    print("number of anomalous elements: ", len(anomalies))
                    print("cutoff percentage: ", format(percentage*100., '.2e'), "%")
            #apply two-tail
            else:
                head = np.max(anomalies[anomalies<(len(x0)/2)])
                tail = np.min(anomalies[anomalies>(len(x0)/2)])
                x = x0[head+1:tail]
                y = y0[head+1:tail]
                percentage = 1 - ((tail - head)/len(x0))
                print("two-tail sanitisation applied.")
                print("number of sanitised elements: ", int(percentage*len(x0)))
                print("number of anomalous elements: ", len(anomalies))
                print("cutoff percentage: ", format(percentage*100., '.2e'), "%")
                
        if percentage > 0.3:
            warnings.warn("significant cutoff is applied. take care of the data quality.", RuntimeWarning)
            
        return x, y
        
    def _remove_offset(self, x:np.ndarray, y:np.ndarray, method:str) -> np.ndarray:
        if method == "mean":
            y -= np.mean(y)
            return x, y
            
        elif method == "filter":
        
            y = signal.sosfilt(signal.butter(2, 1, 'hp', fs=self.sampling_frequency, output='sos'), y)
            
            percent = 0.05
            y = y[int(len(y)*(percent)):int(len(y)*(1-percent))]
            x = x[int(len(x)*(percent)):int(len(x)*(1-percent))]
            x -= min(x)
            
            return x, y
        
        else:
            raise ValueError("only 'mean' and 'filter' methods supported.")
            
    
    def _interpolate(self, x:np.ndarray, y:np.ndarray) -> interpolate.CubicSpline:
        return interpolate.CubicSpline(x, y)
        
        
        
    def _get_spectrum(self, x:np.ndarray, interpolator:interpolate.CubicSpline, density:int) -> np.ndarray:
        x_sample = np.linspace(x[0], x[-1], density)
        yf = np.fft.fftshift(np.fft.fft(interpolator(x_sample)))
        xf = np.fft.fftfreq(len(x_sample))
        xf = np.fft.fftshift(xf)
        yf = np.abs(yf[xf>0])
        xf = xf[xf>0]
        lambda_ = np.abs(np.mean(np.diff(x_sample))/xf)
        
        return lambda_, yf
        
    def calibrate(self, metres_per_microstep:float, density:int, remove_offset_method:str):
        self.metres_per_microstep = metres_per_microstep
        self.x *= self.metres_per_microstep
        
        self.x, self.y = self._auto_sanitise(self.x, self.y)
        self.x, self.y = self._remove_offset(self.x, self.y, remove_offset_method)
        interpolator = self._interpolate(self.x, self.y)
        self.lambda_, self.yf = self._get_spectrum(self.x, interpolator, density)
        
        
        
    def visualise(self, save:bool=False, interactive:bool=False, manual_lim:tuple=None):
        if self.metres_per_microstep is None:
            raise ValueError("calibrate before attempt to visualise.")
        
        if interactive is False:
            fig, ax = plt.subplots(figsize=(9, 6))
            ax.plot(self.lambda_, self.yf)
            ax.set_xlabel("Wavelength (m)")
            ax.set_ylabel("Amplitude")
            
                
            if manual_lim is not None:
                ax.set_xlim(*manual_lim)
                
            
            if save is True:
                fig.savefig('./spectrum_visualisation.png')
                
            return fig, ax
        
        elif interactive is True:
            if manual_lim is not None:
                warnings.warn("x-limiting is not available for interactive visualisation."
                              , RuntimeWarning)
            
            
            df = pd.DataFrame({
                    'x':self.lambda_,
                    'y': self.yf
                })
            label = {'x':'Wavelength (m)', 'y':'Amplitude'}
            fig = px.line(df, x='x', y='y', labels=label)
            if save is True:
                warnings.warn("interactive figures need extra packs to be saved, refer to documentation."+
                              " or you can manually save the static image."
                              , RuntimeWarning)
            fig.show()
            
        else: raise ValueError('only booleans are supported for "interactive" argument.')
        
        
            
            
        
        

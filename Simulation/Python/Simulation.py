import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import warnings
from IPython.display import display
from typing import Union


class InterferogramSimulator():
    def __init__(self, x:np.ndarray, y:np.ndarray, n_mode:int, n_sigma:int):
        self.N_SIGMA = float(n_sigma)
        self.N_MODE = int(n_mode)
        self.X = x
        self.y = y
        self.source_list = []
        
        self._set_plot_style()
        self._set_warning()
        
        if self.X.shape != self.y.shape:
            raise ValueError('mesh shape must be the same.')
    
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
        
    @staticmethod
    def auto_mesh(sampling_frequency:Union[int, float]=50, motor_speed:Union[int, float]=30000
                    , start_position:float=-1e7, conversion_factor:float=1e-11) -> np.ndarray:
        
        end_position = -start_position
        sampling_distance = motor_speed/sampling_frequency
        n_sample = int((end_position - start_position)/sampling_distance)

        x = np.linspace(start_position, end_position, n_sample)*2.0*conversion_factor
        y = np.zeros(len(x))

        return x, y

        
    def _gaussian_spectrum(self, intensity:float) -> np.ndarray:
        spectrum = np.linspace(-self.N_SIGMA, self.N_SIGMA, self.N_MODE)
        spectrum = np.exp(-spectrum*spectrum/2.)
        spectrum /= np.sum(spectrum)
        spectrum *= intensity
        return spectrum
    
    def _adjusted_gaussian_spectrum(self, intensity:float, sigma:float) -> np.ndarray:
        spectrum = np.linspace(-self.N_SIGMA*sigma, self.N_SIGMA*sigma, self.N_MODE)
        spectrum = np.exp(-spectrum*spectrum/2)
        spectrum /= np.sum(spectrum)
        spectrum *= intensity
        return spectrum
        
    def add_gaussian(self, wavelength:float, sigma:float, intensity:float) -> np.ndarray:
        wavelengths = np.linspace(-self.N_SIGMA*sigma, self.N_SIGMA*sigma, self.N_MODE) + wavelength
        spectrum = self._gaussian_spectrum(intensity)
        #spectrum = self._adjusted_gaussian_spectrum(intensity, sigma)
        
        spatial = np.zeros(len(self.X))
        for i in range(len(wavelengths)):
            spatial += spectrum[i]*np.sin(2.*np.pi*self.X/wavelengths[i]) + spectrum[i]
            
        self.y += spatial
        
        self.source_list.append(['gaussian', wavelength, self.N_SIGMA*sigma, intensity])

        return spatial
    
    def add_square(self, start:float, width:float, intensity:float) -> np.ndarray:
        wavelengths = np.linspace(start, start+width, self.N_MODE)
        spectrum = intensity / self.N_MODE
        
        spatial = np.zeros(len(self.X))
        for i in range(len(wavelengths)):
            spatial += spectrum*np.sin(2.*np.pi*self.X/wavelengths[i]) + spectrum
            
        self.y += spatial
        
        self.source_list.append(['square', start+0.5*width, width/2., intensity])

        return spatial
    
    def inspect_simulation(self) -> pd.DataFrame:
        
        if self.source_list is None:
            raise ValueError('source information is not available for C++ simulations.')
            
        else:
            display(pd.DataFrame(self.source_list
                                , columns=['Type', 'Central Wavelength', 'Effective Width', 'Intensity']))

            return pd.DataFrame(self.source_list
                                , columns=['Type', 'Central Wavelength', 'Effective Width', 'Intensity'])

    
    def visualise(self, save:bool=False, interactive:bool=True):
        if interactive is False:
            fig, ax = plt.subplots(figsize=(9, 6))
            ax.plot(self.X, self.y)
            ax.set_xlabel('Distance from null point in meters')
            ax.set_ylabel('Amplitude')

            if save is True:
                fig.savefig('./interferogram_visualisation.png')
                
            fig.show()

            return fig, ax
        
        elif interactive is True:
            df = pd.DataFrame({
                    'x':self.X,
                    'y':self.y
                })
            label = {'x':'Distance from null point in meters', 'y':'Amplitude'}
            fig = px.scatter(df, x='x', y='y', labels=label)
            fig.update_traces(marker=dict(
                size=12,
                symbol='circle'
            ))
            if save is True:
                warnings.warn("interactive figures need extra packs to be saved, refer to documentation."+
                              " or you can manually save the static image."
                              , RuntimeWarning)
            fig.show()
            
        else: raise ValueError('only booleans are supported for "interactive" argument.')
            
    
    def _auto_lims(self) -> tuple:
        
        source_data = pd.DataFrame(self.source_list
                            , columns=['Type', 'Central Wavelength', 'Effective Width', 'Intensity'])
                                   
        source_range = np.array([(source_data['Central Wavelength'] + source_data['Effective Width']).to_numpy(),
                       (source_data['Central Wavelength'] - source_data['Effective Width']).to_numpy()])
                                   
        source_range = (np.max(np.ndarray.flatten(source_range)), np.min(np.ndarray.flatten(source_range)))
        distance = 2*(source_range[0] - source_range[1])
        
        return (source_range[1]-distance, source_range[0]+distance)
            
    
    def show_spectrum(self, save:bool=False, interactive:bool=False, manual_lim:tuple=None):
        
        yf=np.fft.fft(self.y)
        xf=np.fft.fftfreq(len(self.X))

        xf=np.fft.fftshift(xf)
        yf=np.fft.fftshift(yf)

        xf_positive = xf[xf>0]
        lambda_ = np.mean(np.diff(self.X))/xf_positive
        
        yf = np.abs(yf[xf>0])
        
        if interactive is False:
            fig, ax = plt.subplots(figsize=(9, 6))
            ax.plot(lambda_, yf)
            ax.set_xlabel("Wavelength (m)")
            ax.set_ylabel("Amplitude")
            
            if (self.source_list is None) and (manual_lim is None):
                warnings.warn("automatic limit no longer available. use manual_lim argument to set x-limits."
                              , RuntimeWarning)
                
            elif (self.source_list is None) and (manual_lim is not None):
                ax.set_xlim(*manual_lim)
                
            elif (self.source_list is not None) and (manual_lim is not None):
                warnings.warn("automatic limit is recommended for Python simulation."
                              , RuntimeWarning)
                ax.set_xlim(*manual_lim)
                
            else:
                ax.set_xlim(*self._auto_lims())
            
            if save is True:
                fig.savefig('./spectrum_visualisation.png')
                
            return fig, ax
        
        elif interactive is True:
            if manual_lim is not None:
                warnings.warn("x-limiting is not available for interactive visualisation."
                              , RuntimeWarning)
            
            
            df = pd.DataFrame({
                    'x':lambda_,
                    'y': yf
                })
            label = {'x':'Wavelength (m)', 'y':'Amplitude'}
            fig = px.line(df, x='x', y='y', labels=label)
            if save is True:
                warnings.warn("interactive figures need extra packs to be saved, refer to documentation."+
                              " or you can manually save the static image."
                              , RuntimeWarning)
            fig.show()
            
        else: raise ValueError('only booleans are supported for "interactive" argument.')
            
    def load_cpp_simulation_result(self, filepath:str):
        data = np.loadtxt(filepath, delimiter=',').T
        self.X = data[0]
        self.y = data[1]
        self.source_list = None
        print('simulation result loaded. follow the usual routine to inspect.')

# Interferogram Simulation

---

## Python Version

### Dependence Summary
- `numpy`
- `matplotlib` (static visualisation)
- `plotly.express` (interactive visualisation)
- `pandas` (required for data preprocessing for `plotly`)
- `warnings` (for `RuntimeWarning`), `IPython.display` (for summary tables)
- `typing` (for `Union` type hints)
- `seaborn` (optional, for setting plot style)

### Public Methods

---

```python
InterferogramSimulator(x:np.ndarray, y:np.ndarray, n_mode:int, n_sigma:int)
```
Initialiser to instantiate the class `InterferogramSimulator`.

Input:
- `x, y`: np.ndarray. mesh points for the simulation. using the helper function `auto_mesh()` is encouraged.
- `n_mode`: int. number of discretised points used to simulate the continuous spectrum. affect accuracy.
- `n_sigma`: int. number of sigmas taken into account in statics of the simulation. affect accuracy.

Return:
A instance of `InterferogramSimulator`.



```python
@staticmethod
InterferogramSimulator.auto_mesh(sampling_frequency:Union[int, float]=50, motor_speed:Union[int, float]=30000, start_position:float=-1e7, conversion_factor:float=1e-11) -> np.ndarray:
```

Static helper function. Can be used to generate grid.

Input:
- `sampling_frequency`: int/float. sampling frequency of the detector.
- `motor_speed`: int/float/ motor speed.
- `start_position`: float. starting position of the sampling.
- `conversion_factor`: float. the conversion factor from microstep to meter. 

Return:
- `x, y`: np.ndarray. mesh points for simulation. 


```python
InterferogramSimulator.add_gaussian(wavelength:float, sigma:float, intensity:float) -> np.ndarray
```
Adding Gaussian source to the simulation profile.

Input:
- wavelength: float. The central wavelength of the gaussian source.
- sigma: float. standard deviation of the gaussian source.
- intensity: float. relativive intensity in arbitrary unit.

Return:
- spatial: np.ndarray. intensity according to input position information.





```python
InterferogramSimulator.add_square(start:float, width:float, intensity:float) -> np.ndarray
```
Adding square source to the simulation profile.

Input:
- start: float. The starting position of the square profile.
- width: float. width of the square profile.
- intensity: float. relativive intensity in arbitrary unit.

Return:
- spatial: np.ndarray. intensity according to input position information.
    
    


```python
InterferogramSimulator.inspect_simulation() -> pd.DataFrame
```

Inspect the summary of source simulated.

Input: N/A

Return: 
- pd.Dataframe. Contains a dataframe of the simulated source and spectrum information.




```python
InterferogramSimulator.visualise(save:bool=False, interactive:bool=True)
```
Visualise the spatial distribution. 

Input:
- save: boolean. save the produced figure or not. this argument is ignored when `interactive=True`. interactive plots can be saved manually as `.png` file. (to save as `.html`, you will need `kaleido`.)

- interactive: boolean. use `matplotlib` if `False`, `plotly` if `True`.

Return: ignored.




```python
InterferogramSimulator.show_spectrum(save:bool=False, interactive:bool=False, manual_lim:tuple=None)
```

Show the wavelength spectrum reconstructed from Fourier transform.

Input:
- save: boolean. save the produced figure or not. this argument is ignored when `interactive=True`. interactive plots can be saved manually as `.png` file. 

- interactive: boolean. use `matplotlib` if `False`, `plotly` if `True`. 

- manual_lim: tuple. with shape `(lower_limit, upper_limit)`. this argument manually adjust the x-limits of static plot by `matplotlib`. you are advised to use automatic limit by leaving this argument as `None` when the data is not imported from outside.




```python
InterferogramSimulator.load_cpp_simulation_result(filepath:str)
```

Load C++ simualtion results file for visualisation. 

Input:
- filepath: str. this should be the direct output from C++ simulation. `inspect_simulation()` method will be disabled once this method is invoked.

Return: N/A
    
---


## C++ Version

### Dependence Summary
- `eigen-3.4.0` (for array operations)



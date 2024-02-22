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

### Public Methods

- `add_gaussian()`

```python
InterferogramSimulator.add_gaussian(wavelength:float, sigma:float, intensity:float) -> np.ndarray
```
Adding Gaussian source to the simulation profile.

Input:
wavelength: float. The central wavelength of the gaussian source.
sigma: float. standard deviation of the gaussian source.
intensity: float. relativive intensity in arbitrary unit.



## C++ Version

- `eigen-3.4.0` (for array operations)



# Imperial College Year 2 Interferometry Lab Coding

---

There are five major components that are available at the moment. Please see `README.md` in each folder for further instructions. 

All classes take modular design such that all related functionalities are encapsulated in individual methods. Please see the DocStrings in each Python file for reference.

All implemented C++ classes have similar structures compared to Python implementations.

---

See the `Example.py` files for a quick demonstration of usage. C++ examples are included in the `main.cpp` files in the `main` functions. Note the C++ files need to be compiled locally with `cmake`.

---

Using `Example.py`s:
- open terminal, then `cd /path/to/the/python/file`.
- execute by `python Example.py test_data.txt` (remove `test_data.txt` for the `Simulation` file).
- your advised to use these classes in `jupyter notebook` to enjoy the interactive plottings.

---

`quick_plot.py`:
- Content: the quick visualisor for lab data, automatically check which machine has the valid data recorded.
- C++ available?: No.

`Simulations`:
- Content: Simulate interferogram with a very simple interface and interactive visualisor. All methods are optimised by vectorisation.
- C++ available?: Yes. Depends on `Eigen`.

`CrossingPoints`:
- Content: Find the crossing points of interferogram, support butterworth filtering for removing the offset.
- C++ available?: Yes. Depends on `Eigen` and `Iir`.

`GlobalCalibration`:
- Content: Apply global step-to-meters conversion.
- C++ available?: No.

`LocalCalibration`:
- Content: Apply local calibration conversion.
- C++ available?: No.



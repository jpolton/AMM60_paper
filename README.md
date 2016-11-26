# AMM60_paper

**25 Nov**
----------
To many copies of delta_diagnose function and methods. Removed `IT/AMM60_tools.py` and created symbolic link::

    ln -s FASTNEt/amm60_data_tools.py IT/.

Fixed `internaltidemap_AMM60_paper.ipynb` and `pycnocline_mod_obs_virtual_mooring.ipynb` to diagnose delta with function calls. This required treatment of 3D model data or 1D mooring data.
`FASTNEt/amm60_data_tools.py` has therefore been changed (it tests whether the "profile" data is 1D or 3D in space). In particular a **new output field** is added which will break stuff. My stuff is all fixed to work with this.


---
**26 Nov**

* Check `pycnocline_mod_obs_virtual_mooring.ipynb` works
* start edit `diagnose_delta from amm60_data_tools.py`, for running mean
* Implement running mean processing in `diagnose_delta`
* Implement data clipping in `pycnocline_mod_obs_virtual_mooring.ipynb` for pycnocline diagnostics that were computed outside the range that the Doodson filter gave values for no-tide. Note, however, since the running window is computer over +/- 1.5 day window, the first 1.5 days from the first delta_nt will not have from 1.5 increasing to the 3 days of data in the analysis.

---
**To-do**

`internaltidemap_AMM60_paper.ipynb`:
Update to do running mean averaging rather than 3 day chunking
Make some manuscript plots

`pycnocline_mod_obs_virtual_mooring.ipynb`:
Export to py file
Make some manuscript plots

`diagnose_delta from amm60_data_tools.py`:
Karen: make sure it still works for your FFT code

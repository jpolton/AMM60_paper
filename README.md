# AMM60_paper

Removed IT/AMM60_tools.py and created symbolic link
ln -s FASTNEt/amm60_data_tools.py IT/.

Fixed internaltidemap_AMM60_paper.ipynb and pycnocline\_mod\_obs\_virtual\_mooring.ipynb to diagnose delta with function calls. This required treatment of 3D model data or 1D mooring data.
amm60_data_tools.py has therefore been changed. In particular a new output field is added which will break stuff. My stuff is all fixed to work with this.

**To-do**
---------
internaltidemap_AMM60_paper.ipynb:
Update to do running mean averaging rather than 3 day chunking

pycnocline\_mod\_obs\_virtual\_mooring.ipynb:
Update to do running mean averaging rather than 3 day chunking
Export to py file

diagnose\_delta:
Update to do running mean averaging rather than 3 day chunking


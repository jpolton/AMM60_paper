I removed amm60_data_tools.py from this directory.
It should ONLY appear at ../FASTNEt/amm60_data_tools.py

This might break some code in this directory since the path will not work by default. 

Solution- fix path:

	import sys
	sys.path.append('../FASTNEt/') # Add the directory with the amm60_data_tools.py file to path
	import amm60_data_tools # The new AMM60_tools



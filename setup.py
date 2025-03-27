from setuptools import setup, find_packages
from XRBID import __version__

setup(
	name='XRBID', 
	version=__version__,
	url='https://github.com/QosmicQi/XRBID',
	author='Qiana Hunt',
	author_email='qiana.hunt@uleth.ca',
	packages=find_packages(where="XRBID"),
    	package_dir={"": "XRBID"}
	#py_modules=['Align', 'AutoPhots', 'CMDs', 'DataFrameMod', 'Headers', 'HRs', 'ImageSearch', 'Sources', 'WriteScript']
)

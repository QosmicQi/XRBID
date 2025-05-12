(chap:errors)=
# Errors and How to Fix Them

Here is a list of errors I've encountered while running these programs and what I did to fix them. This is by no means a comprehensive list, but you may find it useful all the same.

* When trying to run `import astroquer.mast as Observation`, I encountered an error stating that `StructuredUnit` could not be imported from `astropy.units`.

  - **SOLUTION:** In the command line, run `conda update astropy`. Then, reopen `python` or the `.ipynb` and run again the import.
* While importing `TweakReg` with `from drizzlepack import tweakreg`, particularly on a MAC with a Silicon chip, you may encounter an importation error. This is caused by some problem during the installation process. 

  - **SOLUTION:** In the terminal, try to force conda to reinstall `scipy`, `libgfortran`, and `libgfortran5` with `conda install -c conda-forge --force-reinstall scipy libgfortran libgfortran5`. If this doesn't work, you can try uninstalling and reinstalling `drizzlepac` with `pip uninstall drizzlepac` or `conda remove drizzlepac`, followed by `conda install -c conda-forge -c astroconda drizzlepac`. If that fails, consult with ChatGPT and try its solutions until something works.
* While running `TweakReg`, I encountered an error that read: `OSError: Empty or corrupt FITS file`.

  - **SOLUTION:** At least one of the `FITS` files in your input list is bad. If you run each individual file through `TweakReg`, you may eventually find one (or more) that returns 
    `ValueError: Input file '<filename>' is neither a GEIS file nor a FITS file.` Remove these from your item list. You can also attempt to open each FITS file in a for loop using a try and except block to find the bad files; this is a quicker method.
* When trying to run a bash script obtained by running WriteDS9, received the error `application-specific initialization failed: couldn't connect to display ":0"`. 

  - **SOLUTION:** Try opening an XQuartz terminal either through the application or by typing `xterm &` into the terminal. The bash script should be executable through that new `xterm` terminal. If you receive an additional error `application-specific initialization failed: couldn't connect to display ":1"`, run `export DISPLAY=:0` and try again. One may also try to following steps in the standard command line terminal (and if all else fails, feed your errors into ChatGPT and follow the instructions given):
```
pkill Xquartz
open -a XQuartz
echo $DISPLAY
export DISPLAY=:0
echo 'export DISPLAY=:0' >> ~/.zshrc
source ~/.zshrc
xhost + localhost
xterm &
``` 

* When running the `DS9` bash script from `WriteDS9`, you may encounter an error telling you that `command ds9 not found`, or something to that effect. This happens if `DS9` is not installed as expected. 

  - **SOLUTION:** in the `.sh` script, replace the `ds9` command with the path to your `DS9` installation (e.g. `/Applications/<some path>/ds9`). 
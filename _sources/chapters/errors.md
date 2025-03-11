(chap:errors)=
# Errors and How to Fix Them

Here is a list of errors I've encountered while running these programs and what I did to fix them. This is by no means a comprehensive list, but you may find it useful all the same.

* When trying to run `import astroquer.mast as Observation`, I encountered an error stating that `StructuredUnit` could not be imported from `astropy.units`.
    **SOLUTION:** In the command line, run `conda update astropy`. Then, reopen `python` or the `.ipynb` and run again the import.
* While running `TweakReg`, I encountered an error that read: `OSError: Empty or corrupt FITS file`.
    **SOLUTION:** At least one of the `FITS` files in your input list is bad. If you run each individual file through `TweakReg`, you may eventually find one (or more) that returns 
    `ValueError: Input file '<filename>' is neither a GEIS file nor a FITS file.` Remove these from your item list. You can also attempt to open each FITS file in a for loop using a try and except block to find the bad files; this is a quicker method.
* When trying to run a bash script obtained by running WriteDS9, received the error `application-specific initialization failed: couldn't connect to display ":0"`. 
     **SOLUTION:** Try opening an XQuartz terminal either through the application or by typing `xterm &` into the terminal. The bash script should be executable through that new terminal. If you receive an additional error `application-specific initialization failed: couldn't connect to display ":1"`, run `export DISPLAY=:0` and try again. If all else fails, feed your errors into ChatGPT and follow the instructions given. 
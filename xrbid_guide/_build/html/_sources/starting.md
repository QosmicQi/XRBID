(chap:starting)=
# Getting Started

There is a list of different software you'll need to properly analyze *CXO* and *HST*, which I've compiled here. 

```{note}
If you're using a Mac and your operating system is older than version 12.0, you may encounter unexpected issues with these applications. Please ensure your OS is up to date, if possible, before continuing! This chapter assumes you're using a Mac, but if not, most of these programs will have a Linux installation as well.
```

You will also need a working knowledge of basic `UNIX` to understand some of these instructions. You can find a [helpful cheatsheet of commands here](https://mally.stanford.edu/~sr/computing/basic-unix.html).


```{note}
If you encounter any errors while running these programs, please check Chapter \ref{chap:errors}, **Errors and How to Fix Them**. There I've listed the errors I ran into while creating this document and how I fixed them, so there's likely information in there you may find useful for solving your own issue. If your issue isn't there, it's probably on StackExchange, or you can ask ChatGPT for recommendations.
```

## Installing Software

### `Anaconda` and `python`
The code I use in my work is based on `python`, so you will need to set that up before starting with the data[^1]. The best way to do this is to start with Anaconda. This is generally a painless process, as Anaconda takes care of the rest of your `python` setup for you, except for a few specialized packages you may need.  

If you haven't already, install the appropriate version of Anaconda from here:
https://docs.anaconda.com/anaconda/install/

You'll want to install a command line version rather than the graphic interface (the second tab under the `macOS/Linux installation`). Make sure you're installing the right version for your chip; you can check under `About This Mac` whether your computer has an Intel or an Apple chip. 

Once you install Anaconda, you can use the `conda` command to handle future command line installations and other tasks (such as updating a package). You may also want to use `conda` to set up a specialized environment under which you can work, such as `stenv` (recommended).
`stenv` is a pre-made astronomy-based environment that will automatically install a ton of astronomy-related packages. You can find instructions on setting up `stenv` here:
https://stenv.readthedocs.io/en/latest/getting_started.html

Skip the `conda` installation instructions and follow the instructions under `Choose an stenv release`, making sure to follow the `conda` instructions and not the `micromamba` or `mamba` instructions. To enter the `stenv` environment after setting it up, you'll simply type into the command line (AKA the terminal): 
```
conda activate stenv
``` 

You will see the environment in your terminal change from `(base)` to `(stenv)`. 

Also make sure the `python` distribution that Anaconda installed is up to date (version 3.0 or higher is preferred). You can do this by entering into the terminal:
```
conda update python
```

If you encounter other issues with `conda`, try running: 
```
conda update conda
```

### The Chandra Source Catalog and related tools

To find and download *CXO* sources, you will need the Chandra Source Catalog (`CSC`) desktop app (unless you're using their webpage). This requires `Java`: http://www.java.com

Next, install the `CSC` desktop app, `CSCview`: http://cda.cfa.harvard.edu/cscview/

`CSC` outputs data in the form of `VOTables` or `TSV` (tab-separated value) tables, which are basically forms of spreadsheets. You may find it easier to open these formats using a program like `TOPCAT`: https://www.star.bris.ac.uk/~mbt/topcat/#install 

I recommend running the `Standalone Jar File`, but you can also do a full install of the graphic interface with: 
```
url -OL http://www.starlink.ac.uk/topcat/topcat-all.dmg
```

### Imaging with `SAOImageDS9`

It is also incredibly useful to have some visualization tool on hand that can display color images of the galaxy and plot sources from coordinates. I use `DS9`, although there is a new version of [CARTA](https://cartavis.org/) that is going to be released soon that will have the same capabilities as `DS9` in a more functional, cleaner program. Until then, `DS9` is the way to go: https://sites.google.com/cfa.harvard.edu/saoimageds9

Unless you know what you're doing, I recommend installing the Aqua version for Mac, which will allow you to open `DS9` like any other app. This is found by clicking the large `Download` button, selecting the MacOS drop-down menu, and reading through the instructions[^2]. Note, if the Aqua version gives you a "damaged application" error when you try to open `DS9`, you will also have to run the following code to fix it: 

```
xattr -c /Applications/SAOImageDS9.app
``` 

If you decide to install the command line version, then you will need to follow the instructions for `Darwin X11`. The home directory into which you want to move the `DS9` files will probably look like: ```/usr/local/bin```. Navigate to this directory to ensure it's there. If your `/usr/local/` directory is empty, add the `bin` directory with: 
```
mkdir bin
```

If while installing `DS9` you encounter an error that states the program has an unknown developer, you will need to change the security settings of your computer to allow you to proceed anyway. If you have no issues installing `DS9` but have trouble opening it, then make sure `Darwin X11` and `XQuartz` are installed and up to date on your computer.

### Other `python` packages

If you intend to query *HST* images (or any publicly-available archival astronomy data) through `python`, you will likely need to install `PyVO`, which may be found here: 
https://pyvo.readthedocs.io/en/latest/  

If you plan on creating mosaics from the *HST* images, you will also want to install `astroquery`, found here: https://astroquery.readthedocs.io/en/latest/#using-astroquery

To make mosaics, you'll also need to install `DrizzlePac`[^3], which contains the `AstroDrizzle` function. *But* `DrizzlePac` will probably give you issues if installed on its own, so it's recommended to install and run it in a custom environment (`stenv`). `stenv` should have automatically installed `DrizzlePac`, but in the event it didn't, you can find it here: 
https://drizzlepac.readthedocs.io/en/latest/

Remember to enter the `stenv` environment first before running the installation: 
```
conda activate stenv
```  

Some of these packages may require `git` to install. Most Mac computers will have this software already installed, but if you're using an older machine, you may have to install it manually: 
https://github.com/git-guides/install-git

## Starting `python`

I've found it easiest to run the analysis out of a `iPython` notebook, which allows you to edit and run individual cells of code so that you can repeat each step of the analysis independently. To open `python` in as a notebook, open your command line and enter: 
```
jupyter notebook &
```
(where the optional `&` opens the notebook separately from the command line, keeping the command line free for continued use). If you wish to run it out of the `stenv` environment, be sure to run `conda activate stenv` first. 

This command will open an interface in a new browser window/tab with a list of files in the directory from withing which you called the command. You can select a pre-existing notebook (which have the suffix `.ipynb`), or create a new one using the `New > Notebook` option in the menu on the right. This will open a new `.ipynb` for you to work within. 

Alternatively, you can run this entire project out of Google Colaboratory, a `python` compiler that runs out of Google Drive. This has the benefit of being able to run from anywhere from any computer, regardless of the operating system. However, the downside is that it's much slower than running `python` on your laptop and is prone to timing out if too much memory is used at once. To open `python` in Google, enter your Google Drive, navigate to your preferred directory, select `+ New > More > + Connect more apps` and find Google Colab on the list (it has an orange CO logo). Afer installing it once, you'll be able to find the app under `+ New > More` without having to install it again. From there, you can install all necessary astronomy packages into Google Colab using `!pip install` within one of the cells. You should only have to do this once per package, which should then be accessible to any other Google Colab notebook you make. 

## Starting `CSC` (and complementary programs)

To open `CSC` on a Mac, you should be able to click on the application icon and have it open directly. Otherwise, you may call it from the command line from the directory where it was installed:  

```
java -jar cscview.jar
``` 
or, if you have a newer version of `java` than I use: 

```
java --add-modules java.se.ee -jar cscview.jar
``` 

Likewise, to view the `CSC` data table in a user-friendly format, you can open `TOPCAT` by calling it in the command line from the directory where it was installed:  
```
java -jar topcat-*.jar  # Whatever your version of topcat is saved as
```

From there, you should be able to navigate to your desired table and open it to view. 

## A Note About Data Management

Throughout this guide, I will be organizing my data as Pandas `DataFrames`. I've found this is the easiest way to view large collections of data. There are certainly more efficient data formats, but I appreciate the ability to seamlessly mix strings and numbers into a readable (and writeable) table. One thing to note, however: `DataFrames` are not very good at saving numerical data to a high degree of precision, so be wary when reading and writing values. If you notice minor discrepancies in your numbers when working on your data on different days, the `DataFrames` may be the culprit! 

At the end of the guide, I will provide the code I use to manage and manipulate these DataFrame files, which may be of use to you as well or may form the baseline of your own code. These are found in Section \ref{sec:script-dataframemod}.


[^1]: Alternatively, {cite}`chandar20` choose to use `IRAF` for their photometric measurements. I will not go into detail about how that is done, but if you find you'd rather use `IRAF` (or `PyRAF`), [more information about its installation can be found here](https://faculty1.coloradocollege.edu/~sburns/courses/18-19/pc362/Anaconda_IRAF_install.html).
[^2]: If your OS does not support any of these versions of `DS9`, you can find archival versions here: https://ds9.si.edu/archive
[^3]: https://www.stsci.edu/scientific-community/software/drizzlepac
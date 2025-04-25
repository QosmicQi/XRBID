(chap:starting)=
# Getting Started

There is a list of different software you'll need to properly analyze *CXO* and *HST*, which I've compiled here. 

```{important}
If you're using a Mac and your operating system is older than version 12.0, you may encounter unexpected issues with these applications. Please ensure your OS is up to date, if possible, before continuing! This chapter assumes you're using a Mac, but if not, most of these programs will have a Linux installation as well.
```

You will also need a working knowledge of basic `UNIX` to understand some of these instructions. You can find a [helpful cheatsheet of commands here](https://mally.stanford.edu/~sr/computing/basic-unix.html).


```{note}
If you encounter any errors while running these programs, please check {ref}`chap:errors`. There I've listed the errors I ran into while creating this document and how I fixed them, so there's likely information in there you may find useful for solving your own issue. If your issue isn't there, it's probably on StackExchange, or you can ask ChatGPT for recommendations.
```

## Installing Software

### `Anaconda` and `python`
The code I use in my work is based in `python`, so you will need to set that up before starting with the data[^1]. The best way to do this is to start with Anaconda. This is generally a painless process, as Anaconda takes care of the rest of your `python` setup for you, except for a few specialized packages you may need.  

If you haven't already, install the appropriate version of Anaconda from here:
https://docs.anaconda.com/anaconda/install/

You'll want to install a command line version rather than the graphic interface (the second tab under the `macOS/Linux installation`). Make sure you're installing the right version for your chip; you can check under `About This Mac` whether your computer has an Intel or an Apple chip. 

Once you install Anaconda, you can use the `conda` command to handle future command line installations and other tasks (such as updating a package). You may also want to use `conda` to set up a specialized environment under which you can work, such as `stenv` (recommended).
`stenv` is a pre-made astronomy-based environment that will automatically install a ton of astronomy-related packages. You can find instructions on setting up `stenv` here:
https://stenv.readthedocs.io/en/latest/getting_started.html

Follow the instructions under `Choose an stenv release` (including downloading the proper release for your environment [here](https://github.com/spacetelescope/stenv/releases)), making sure to follow the `conda` instructions (third tab) and not the `micromamba` or `mamba` instructions. To enter the `stenv` environment after setting it up, you'll simply type into the command line (AKA the terminal): 
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

### Other `python` packages

#### `XRBID` and `git`
Some of these packages may require `git` to install. Most Mac computers will have this software already installed, but if you're using an older machine, you may have to install it manually: 
https://github.com/git-guides/install-git

By now, you should have already installed the `XRBID` package from GitHub. If you haven't, you can find the installation instructions at this repository: https://github.com/QosmicQi/XRBID

You can either install this in the `stenv` environment, or make a new one to make sure this package doesn't accidentally interfere with your normal `python` setup. For example, you can make a copy of the `stenv` environment and install XRBID there: 
```
conda create --name xrbenv --clone stenv
conda activate xrbenv
```

```{note}
I am currently making changes to the installation and read-in of the modules in XRBID, so some of the code in the rest of this guide may soon change. Keep an eye out for XRBID version 2.0!! 
```

#### `PyVO`
If you intend to query *HST* images (or any publicly-available archival astronomy data) through `python`, you will likely need to install `PyVO`, which may be found here: 
https://pyvo.readthedocs.io/en/latest/  

#### `astroquery`
If you plan on creating mosaics from the *HST* images, you will also want to install `astroquery`, found here: https://astroquery.readthedocs.io/en/latest/#using-astroquery

#### `DrizzlePac`
To make mosaics, you'll also need to install `DrizzlePac`[^3], which contains the `AstroDrizzle` function. *But* `DrizzlePac` will probably give you issues if installed on its own, so it's recommended to install and run it in a custom environment (`stenv` or `xrbenv`). `stenv` should have automatically installed `DrizzlePac`, but in the event it didn't, you can find it here: 
https://drizzlepac.readthedocs.io/en/latest/

Remember to enter the proper environment first before running the installation, for example: 
```
conda activate stenv
```  

### The Chandra Source Catalog and related tools

To find and download *CXO* sources, you will need the Chandra Source Catalog (`CSC`) desktop app (unless you're using their webpage). This requires `Java`: http://www.java.com

```{note}
If you're installing the Linux version of Java, you should move the downloaded `.tar` file to a directory called `/usr/java/` and open it there. Then, to create a simlink that can be use anywhere on your computer, run `sudo ln -s /usr/java/jre1.<distribution version>/bin/java /usr/local/bin/java`. This should cause `java` to show up in your `/usr/local/bin` file in the same color as any other simlinks you created. 
```

Next, install the `CSC` desktop app, `CSCview`: http://cda.cfa.harvard.edu/cscview/

#### `TOPCAT`
`CSC` outputs data in the form of `VOTables` or `TSV` (tab-separated value) tables, which are basically forms of spreadsheets. You may find it easier to open these formats using a program like `TOPCAT` (thought this is not strictly required): https://www.star.bris.ac.uk/~mbt/topcat/#install 

You can either run the `Standalone Jar File` (which will require you to open the program through the terminal), or you can do a full install of the graphic interface by downloading the `PKG Packaged Application` (currently only available on Mac OS X) or by entering into the terminal: 
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

```{note}
If installing DS9 or other apps on a Linux, you can create a symlink to open the app from the command line by moving the downloaded application to the `/usr` directory and then creating the link with: 
`sudo ln -s /usr/<app_name> /usr/local/bin/<app_name>`
```

If while installing `DS9` you encounter an error that states the program has an unknown developer, you will need to change the security settings of your computer to allow you to proceed anyway. If you have no issues installing `DS9` but have trouble opening it, then make sure `Darwin X11` and `XQuartz` are installed and up to date on your computer.

Please note, there are some process in this guide that open `DS9` as a bash script (e.i. `XRBID.WriteScript.WriteDS9`). These scripts require you open a special terminal before running them, or you may run into errors with `Darwin X11` or `XQuartz`: 
```
xterm &
```

## Starting `python`

If everything is installed properly, you should be able to start up `python` from the command line by simply typing: 

```
python
```

When it opens, it will print some version and licensing information. The cursor where each new line of code starts will follow after `>>>` instead of `(base)` or `(stenv)`. 

I've found it easiest to run the analysis out of a `iPython` notebook, which allows you to edit and run individual cells of code so that you can repeat each step of the analysis independently. To open `python` in as a notebook rather than in the command line, open your command line and enter: 
```
jupyter notebook &
```
(where the optional `&` opens the notebook separately from the command line, keeping the command line free for continued use). 


This command will open an interface in a new browser window/tab with a list of files in the directory from withing which you called the command. You can select a pre-existing notebook (which have the suffix `.ipynb`), or create a new one using the `New > Notebook` option in the menu on the right. This will open a new `.ipynb` for you to work within. Make sure the kernel matches the environment you wish to run the code out of. 


```{note}
If you wish to run it out of the `stenv` environment, be sure to run `conda activate stenv` first, regardless of whether you're running your code out of the command line or an `iPython` notebook. 
```

Alternatively, you can run this entire project out of Google Colaboratory, a `python` compiler that runs out of Google Drive. This has the benefit of being able to run from anywhere from any computer, regardless of the operating system. However, the downside is that it's much slower than running `python` on your laptop and is prone to timing out if too much memory is used at once. To open `python` in Google, enter your Google Drive, navigate to your preferred directory, select `+ New > More > + Connect more apps` and find Google Colab on the list (it has an orange CO logo). Afer installing it once, you'll be able to find the app under `+ New > More` without having to install it again. From there, you can install all necessary astronomy packages into Google Colab using `!pip install` within one of the cells.

```{note}
Google Colab does not currently support persistent installations, meaning you will need to re-install all special packages each time you start a new runtime. This includes your installation of `XRBID`!
``` 

## Starting `CSC` (and complementary programs)

To open `CSC` on a Mac, you should be able to click on the application icon and have it open directly. If you chose to download the `Standalone Jar File`, you may call it from the command line from the directory where it was installed:  

```
java -jar cscview.jar
``` 
or, if you have a newer version of `java` than I use, you can call it with: 

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

Provided in the `XRBID` package is a custom module for searching and modifying the `DataFrames` created in this guide called `DataFrameMod()`. The most important tool within this module is the `Find()` function (`XRBID.DataFrameMod.Find`), which allows the user to find entries with specific values under a given header. For example, one may use: 

```
from XRBID.DataFrameMod import Find
Find(df_example, "CSC ID = 2CXO J140312.5+542056")
```
to find a specific X-ray source under the header "CSC ID" in a `DataFrame` called "df_example", or 
```
Find(df_example, ["CSC ID = 2CXO J140312.5+542056", "Theta < 3"])
```
to stack search criteria. Some more examples of its use can be found throughout the guide. 


[^1]: Alternatively, {cite}`chandar20` choose to use `IRAF` for their photometric measurements. I will not go into detail about how that is done, but if you find you'd rather use `IRAF` (or `PyRAF`), [more information about its installation can be found here](https://faculty1.coloradocollege.edu/~sburns/courses/18-19/pc362/Anaconda_IRAF_install.html).
[^2]: If your OS does not support any of these versions of `DS9`, you can find archival versions here: https://ds9.si.edu/archive
[^3]: For more information on `drizzlepac`, visit: https://www.stsci.edu/scientific-community/software/drizzlepac
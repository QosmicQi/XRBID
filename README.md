# About XRBID
XRBID is a custom `python` package created primarily to facilitate the identification of *HST* optical sources associated with X-ray binaries observed by *Chandra*. 

For more information on this process or how these packages are used for practical research, visit the [XRBID Guidebook](https://qosmicqi.github.io/XRBID/chapters/intro.html).

# Installing XRBID in a terminal

From the directory where you would like to download XRBID, run in your command line terminal: 

```
git clone https://github.com/QosmicQi/XRBID
```

Navigate into the XRBID directory and install the modules through the setup file: 
```
cd XRBID
sudo python setup.py install  
```

This should install XRBID and all of the modules it contains into python. Now from python, you should be able to call on these modules, i.e.
```
python
from XRBID import Sources
from XRBID.DataFrameMod import Find
# etc.
```

# Installing XRBID on Google Colab

If you are running your analysis through Google Drive with Google Colab, then you can install the necessary packages by opening a new Colab `iPython` notebook and running the following in an empty cell: 

```
!git clone https://github.com/QosmicQi/XRBID.git
%cd XRBID
!pip install -e .
```

By default, this will clone the repository to a file on the `content` directory with the path `/content/XRBID`. Test the installation by importing the module and its functions into your notebook. 

**NOTE:** Google Colab does not currently have a method for permanently installing GitHub repos, so you will need to rerun the installation every time you reboot Google Colab. You will need to navigate back to `/content/XRBID` with each new session and run the installation command. 

# Updating XRBID
As this module is under active development, you should update it from time to time, especially following bug fixes. In your terminal, run: 

```
pip install --upgrade git+https://github.com/QosmicQi/XRBID.git
```

If you're using Google Colab, then you will need to delete the current cloned repository and reinstall it as before. 

```
!rm -rf /content/XRBID
!cd /content
!git clone https://github.com/QosmicQi/XRBID.git
%cd XRBID
!pip install -e .

```

If that doesn't work, you can always uninstall XRBID and start from scratch using the instructions above. 

# Citing XRBID
If you happen to use this code for research that results in a paper, please add it to the Software section and/or cite it as follows: 

```
Hunt, Q. (2025). XRBID: A Python package for X-ray binary identification (v1.0) [Software]. GitHub. https://github.com/QosmicQi/XRBID
```

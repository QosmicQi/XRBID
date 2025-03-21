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

Test the installation by importing the module and its functions into your notebook. 

# Updating XRBID
As this module is under active development, you should update it from time to time, especially following bug fixes. In your terminal, run: 

```
pip install --upgrade git+https://github.com/QosmicQi/XRBID.git
```

or, if you're using Google Colab: 

```
!pip install --upgrade git+https://github.com/QosmicQi/XRBID.git
```

If that doesn't work, you can always uninstall XRBID and start from scratch using the instructions above. 

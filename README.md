# About XRBID
XRBID is a custom `python` package created primarily to facilitate the identification of *HST* optical sources associated with X-ray binaries observed by *Chandra*. 

For more information on this process or how these packages are used for practical research, visit the [XRBID Guidebook](https://qosmicqi.github.io/XRBID/chapters/intro.html)

# Installing XRBID in `python`

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
```

etc. 

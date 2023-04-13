# Copyright (C) 2023  Sergio Frasca
#  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
Introduction

SnagPy is a Python package for data analysis mainly in the gravitational signal field. 
It is a quasi-porting of the Snag Matlab toolbox developed by the Virgo Rome Data Analysis group.

To start:

-	Put the files contained in the zip where you want, maintaining the folder structure.
-	Set the PYTHONPATH in the environment variables (how to do this depends on the operative system. 
        See for example How to set python path (net-informations.com) 
        or Using PYTHONPATH — Functional MRI methods (bic-berkeley.github.io),
        or anything similar googling "set path python”).
-	Use another folder for your python work.
-	Set the environment variables SNAGPY_PATH, ANT_TAB, SOUR_TAB and HINJ_TAB to the relative files.

External packages to install:

* NumPy
* SciPy
* Matplotlib
* Astropy, jplephem, Skyfield
* Mat73         for reading mat v7.3 files
* h5py          for supporting HDF5 format
* silx

SnagPy is composed by various modules, divided in main modules and secondary modules.

The main modules are in the folder SnagPy, the secondary modules are in the folders 
* Exper, containing experimental modules, 
* Pers, containing the personal or personally modified modules, 
* Examples, containing analysis examples using SnagPy.

In this guide only the main modules are described, and some examples are reported.

The main modules are:

-> GD           the management of the basic container class for 1-D data
-> GD2          for the 2-D data GD
-> DS           data stream management
-> MGD          multiple GD management
-> BASIC        general management routines
-> SERV         service programming routines
-> ML_PY        for using Matlab objects
-> STAT         statistics
-> SIGNAL       signal processing routines
-> IMAGE        image processing
-> ASTROTIME    time and astronomy routines
-> GWDATA       gravitational wave data (sources, antennas, periods…)
-> PSS          periodic source search 
-> BSD          BSD analysis
-> GWOTH        other gw analysis
-> GUISNAG      guis for SnagPy
-> EXT_PACK     module using critical external packages
-> MAN_SUPER    management and supervision
-> PARGPU       parallel and GPU computing
-> WEB_SNAG     routines for the web and phone and tablets apps
-> DEEPSNAG     deep learning for SnagPy
-> FANCY_FIG    fancy figures and other

The secondary modules are in the folders:

* Projects       containing projects related modules
* Exper          containing experimental and test modules
* Pers           containing personal modules

To start with all the basic import in the python shell,
from the SnagPy folder give at the command prompt:
        python -i starting.py

Other documentation:

> Tutorials                     -> dummy_tut
> Programming Guide             -> dummy_pg
> Documentation development     -> dummy_doc

Some hints:

> starting.py contains a "standard" sequence of start up commands.
    Possible use: Starting from the SnagPy folder, to start Python with the start up commands,
    give at the command prompt 
        python -i starting.py

> snippets.py contains snippets of useful code
'''

def dummy_tut():
    '''
    '''

def dummy_pg():
    '''
    '''

def dummy_doc():
    '''
    '''
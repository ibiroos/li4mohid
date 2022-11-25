# LI4MOHID

LI4MOHID is a free and open source [QGIS](https://qgis.org/) plug-in to help researchers and practitioners with their
[MOHID Lagrangian](http://www.mohid.com/pages/models/mohidlagrangian/mohid_lagrangian_home.shtml) modelling workflow. 
LI4MOHID can be used to  input data, run simulations, and visualize results. 

This plug-in has been developed under [MyCOAST](http://www.mycoast-project.org) project.

<p align="center"><img src="./doc/assets/li4mohid_SC.png"></p>
<p align="center">A screenshot of LI4MOHID in QGIS</p>

## Requierements

In order to use the plug-in, is required to have installed:

* **QGIS:** QGIS 3.22 LTR or up is recommended.
* **Lagrangian MOHID:** The lagrangian MOHID model to be executed by the plug-in. The last release
can be downloaded from [MOHID Lagrangian Github site](https://github.com/Mohid-Water-Modelling-System/MOHID-Lagrangian/tags).
  To install it, just download the release folder to your local computer.
  
## Installation


1. Install VTK library for Python QGIS. To install them we recommend the following steps:
    1) Open OSGeo4W Shell
    2) Use next command: pip install VTK
       
  *Note: NCDF4 library is now in the OSGeo4W python package.*

2. Install the plug-in.

Download the li4mohid folder  and copy to C:\OSGeo4W64\apps\qgis\python\plugins folder. 
Open QGIS, click tab “Complements/Manage and install plugins”, select “LI4MOHID” and activate.
An “oil-spill” icon must appear in the toolbar.

## Usage

For usage, please, consult the main [documentation](./doc/MyCoast_LI4MOHID_EN.pdf) 


## History

*2022/11/15:* It was updated to QGIS 3.28 and python 3.9. Use with Lagrangian MOHID v 2022 was checked

## Credits
* *author:*
  + Carlos F. Balseiro (4Gotas, cfbalseiro@4gotas.com)
  + Pedro Montero (INTECMAR, pmontero@intecmar.gal)
* *license:* Copyright (c) INTECMAR 2020. Licensed under MIT
* *funding:* MYCOAST  Interreg Atlantic Programme, Project nr. EAPA 285/2016
             http://www.mycoast-project.org
  <p align="center"><img src="./doc/assets/MyCoast_Logo_small.png"></p>

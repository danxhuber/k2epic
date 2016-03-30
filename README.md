# k2epic
Generate an Ecliptic Plane Input Catalog (EPIC) for the K2 Mission

###Documentation:
- https://archive.stsci.edu/k2/epic.pdf <br/> 
- Huber et al. (2016, ApJS, in press; http://arxiv.org/abs/1512.02643)


###Dependencies:
* Python: numpy,astropy,astroquery <br/> 
* IDL/GDL: ASTROLIB, Coyote Library <br/> 
* Code is tested and functional with GDL 0.9.5 (open source implementation of IDL)

###Usage:

Calling sequence from IDL/GDL:
```
k2epic,ra=inputRA,dec=inputDEC,inrad=inputrad,step=step,epicid=inputepicno
```

**Required input parameters:**	<br/>
inputRA  ... input right ascension (decimal degrees) <br/>
inputDEC ... input declination (decimal degrees) <br/>
inputrad ... input radius (decimal degrees) <br/>
	
**Optional input parameters:**	<br/>
inputepicno	... first EPIC ID for catalog (default=0) <br/>
step 		    ... stepsize in degrees when downloading catalogs from Vizier (default=1) [if the python script times out when download catalogs from Vizier, try to adjust this number]

**Output:** <br/>
Output is stored in the output/ directory, and filenames start with raXXX_deXXX_rXXX, where XXX are the provided input coordinates and radius. Files stored in output/ are: <br/>

* raXXX_deXXX_rXXX.dmc.dat	...	ascii output in DMC (MAST) delivery format
* raXXX_deXXX_rXXX.soc.dat	...	ascii output in SOC (K2 Science Office) delivery format
* raXXX_deXXX_rXXX.fits		  ...	catalog in ascii fits format
* raXXX_deXXX_rXXX.dmc.dat.readme ...	ascii file describing the column headers in the DMC delivery file. 

**Examples:** <br/>
Make a 9 deg radius EPIC around campaign 4 starting with EPIC ID 210282561:
```	
k2epic,ra=59.0758,dec=18.6606,inrad=9,step=1,epicid=210282561
```

Make a 1 deg radius EPIC centered on campaign 1:
```
k2epic,ra=173.94,dec=1.42,inrad=1,step=1,epicid=0
```

Make a 2 deg radius EPIC centered on campaign 2:
```
k2epic,ra=246.13,dec=-22.44,inrad=1,step=1
```

**Details of main subroutines:**
* k2epic.pro is a wrapper for the three main parts of the code, which can also be run separately (by commenting out relevant lines in k2epic.pro)

* getcatalogs.py: python script to download catalogs from Vizier; if anything breaks during the catalog download from Vizier, this is the place to look

* makecat.pro: IDL/GDL code to cross-match downloaded catalogs and calculate Kepler magnitudes

* makeascii.pro: IDL/GDL code to generate ascii output needed for SOC and DMC delivery


**Additional Test Routines:**
	
* testepic.pro: performs basic sanity checks of the produced catalogs, including comparisons to POSS images. Requires ds9 to be installed. See file header fordetails.
	
* epicpropermotion.pro: code to check high proper motion stars in a K2 target list; See file header for details.
	

###Checkk2fov:

Contains a hack of Tom Barclay's k2fov to overplot an EPIC catalog with a K2 FOV. This check is standard procedure to ensure that a generated EPIC catalog overlaps with the campaign coordinates.

To run: <br/> 
1) Edit input.txt and specify some coordinates near the FOV for a given campaign <br/> 
2) run 
```
python checkfov.py input.txt campaignnumber epicpath 
```
where <br/>
campaignnumber = integer specifying campaign number <br/>
epicpath = full path to EPIC catalog in DMC (MAST) format (i.e. a *dmc.mrg file) <br/>

Example using the catalog generated above:
```
python checkfov.py input.txt 1 ../output/ra173.9_de1.4_r1.0.dmc.dat
```
	
The output will be a targets_fov.png plot which shows the overlap between the EPIC and the K2 FOV. Only every 100th source in the catalog is plotted to speed things up. Note that the coordinates future campaigns will have to be edited manually in K2onSilicon.py. These should be the same coordinates that Tom uses for the K2fov code in his github repo.

############################################################################
##### python script to download stellar catalogs from Vizier for K2 fields
############################################################################

###input (arguments passed to script):

##search radius
#rad = 8.5

##stepsize for tiles (all catalogs except Hipparcos + Tycho)
#step = 1.

##coordinates (in degrees)

# campaign 0
#sra = 98.296250
#sde = 21.587778

# campaign 1
#sra = 173.93958
#sde = 1.4172222


###output: binary fits files containing the catalogs for further processing

import numpy as np
import requests
import time
import astropy 
import argparse

import astropy.coordinates as coord
import astropy.units as u

from astropy.io import fits
from astropy.coordinates import Angle
from astroquery.vizier import Vizier

#import matplotlib.pyplot as plt
#import pylab

import pdb

parser = argparse.ArgumentParser(description='Query Vizier to get EPIC input catalogs')
parser.add_argument('sra',help='Input Right Ascension', type=float)
parser.add_argument('sde',help='Input Declination', type=float)
parser.add_argument('rad',help='Search Radius', type=float)
parser.add_argument('step',help='Tile Step Size', type=float)

args = parser.parse_args()
sra=args.sra
sde=args.sde
rad=args.rad
step=args.step

st = (sra-rad)
ntiles=0
while st <= sra+rad:
	st = st + step
	ntiles += 1
	

##########################################################
### Hipparcos
print('Querying Hipparcos')
t0 = time.time()
v = Vizier(columns=['_RAJ2000', '_DEJ2000','HIP', 'Plx', 'e_Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Hpmag', 'e_Hpmag'])

# for some reason we need this first to get more than 50 entries 
Vizier.ROW_LIMIT = -1
v.ROW_LIMIT = -1

result = v.query_region(coord.ICRS(ra=sra, dec=sde, unit=(u.deg, u.deg)), radius=Angle(rad, "deg"), catalog='I/311')

#print(result)
#print(sra,sde,rad)
#pdb.set_trace()

table=result[0]
#table.colnames
#pdb.set_trace()	
#x=np.array(table[:]['Plx'])
#y=np.array(table[:]['e_Plx'])

ra = np.array(table[:]['_RAJ2000'])
dec = np.array(table[:]['_DEJ2000'])
dis = np.sqrt( np.power((ra-sra),2) + np.power((dec-sde),2) )
keep = np.where(dis < rad)

ra = fits.Column(name='_RAJ2000', format='D', array=np.array(table[keep[0]]['_RAJ2000']))
dec = fits.Column(name='_DEJ2000', format='D', array=np.array(table[keep[0]]['_DEJ2000']))
hip = fits.Column(name='HIP', format='D', array=np.array(table[keep[0]]['HIP']))
plx = fits.Column(name='Plx', format='D', array=np.array(table[keep[0]]['Plx']))
e_plx = fits.Column(name='e_Plx', format='D', array=np.array(table[keep[0]]['e_Plx']))
pmra = fits.Column(name='PMRA', format='D', array=np.array(table[keep[0]]['pmRA']))
e_pmra = fits.Column(name='E_PMRA', format='D', array=np.array(table[keep[0]]['e_pmRA']))
pmde = fits.Column(name='PMDE', format='D', array=np.array(table[keep[0]]['pmDE']))
e_pmde = fits.Column(name='E_PMDE', format='D', array=np.array(table[keep[0]]['e_pmDE']))

cols = fits.ColDefs([ra,dec,hip,plx,e_plx,pmra,e_pmra,pmde,e_pmde])
tbhdu = fits.new_table(cols)
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('catalogs/hip.fits',clobber=True)
#print(time.time()-t0)

#plt.plot(x,y,'ro')
#pylab.show()


##########################################################
###### Tycho
print('Querying Tycho')
t0 = time.time()
v = Vizier(columns=['_RAJ2000', '_DEJ2000','TYC1', 'TYC2', 'TYC3', 'pmRA', 'e_pmRA', \
'pmDE', 'e_pmDE', 'BTmag', 'e_BTmag', 'VTmag', 'e_VTmag', 'HIP' ])

Vizier.ROW_LIMIT = -1
v.ROW_LIMIT = -1

result = v.query_region(coord.ICRS(ra=sra, dec=sde, unit=(u.deg, u.deg)), \
radius=Angle(rad, "deg"), catalog='I/259')

table=result[0]
#print(table.colnames)
#pdb.set_trace()	
#x=np.array(table[:]['Plx'])
#y=np.array(table[:]['e_Plx'])

ra = np.array(table[:]['_RAJ2000'])
dec = np.array(table[:]['_DEJ2000'])
dis = np.sqrt( np.power((ra-sra),2) + np.power((dec-sde),2) )
keep = np.where(dis < rad)

ra = fits.Column(name='_RAJ2000', format='D', array=np.array(table[keep[0]]['_RAJ2000']))
dec = fits.Column(name='_DEJ2000', format='D', array=np.array(table[keep[0]]['_DEJ2000']))
tyc1 = fits.Column(name='TYC1', format='20A', array=np.array(table[keep[0]]['TYC1']))
tyc2 = fits.Column(name='TYC2', format='20A', array=np.array(table[keep[0]]['TYC2']))
tyc3 = fits.Column(name='TYC3', format='20A', array=np.array(table[keep[0]]['TYC3']))
pmra = fits.Column(name='PMRA', format='D', array=np.array(table[keep[0]]['pmRA']))
e_pmra = fits.Column(name='E_PMRA', format='D', array=np.array(table[keep[0]]['e_pmRA']))
pmde = fits.Column(name='PMDE', format='D', array=np.array(table[keep[0]]['pmDE']))
e_pmde = fits.Column(name='E_PMDE', format='D', array=np.array(table[keep[0]]['e_pmDE']))
btmag = fits.Column(name='BTMAG', format='D', array=np.array(table[keep[0]]['BTmag']))
e_btmag = fits.Column(name='E_BTMAG', format='D', array=np.array(table[keep[0]]['e_BTmag']))
vtmag = fits.Column(name='VTMAG', format='D', array=np.array(table[keep[0]]['VTmag']))
e_vtmag = fits.Column(name='E_VTMAG', format='D', array=np.array(table[keep[0]]['e_VTmag']))
hip = fits.Column(name='HIP', format='D', array=np.array(table[keep[0]]['HIP']))

cols = fits.ColDefs([ra,dec,tyc1,tyc2,tyc3,pmra,e_pmra,pmde,e_pmde,btmag,e_btmag,\
vtmag,e_vtmag,hip])
tbhdu = fits.new_table(cols)
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('catalogs/tycho.fits',clobber=True)
#print(time.time()-t0)

#pdb.set_trace()	



##########################################################
###### UCAC4

st = (sra-rad)
ix=0

while st <= sra+rad:

	print('Querying UCAC, tile '+str(ix+1)+' of '+str(ntiles))

	# with version 2.3, astroquery apparently lost the ability to correctly query for specific columns
	# instead, we now query for all columns

#	v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'UCAC4', 'of', 'db', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE',\
#	'Tycho-2', '2Mkey', 'Jmag', 'e_Jmag', 'q_Jmag', 'Hmag', 'e_Hmag', 'q_Hmag', 'Kmag',\
#	'e_Kmag', 'q_Kmag', 'Bmag', 'e_Bmag', 'f_Bmag', 'Vmag', 'e_Vmag', 'f_Vmag', 'gmag', \
#	'e_gmag', 'f_gmag', 'rmag', \
#	'e_rmag', 'f_rmag', 'imag', 'e_imag', 'f_imag', 'H', 'LEDA'])

	v = Vizier(columns=['all'])
	#v = Vizier(columns=[':'])
	#v = Vizier(columns=['e_pmRA','e_pmDE','Tycho-2', '2Mkey','e_Jmag', 'q_Jmag', 'e_Hmag', 'q_Hmag', 'e_Kmag', 'q_Kmag', 'e_Bmag', 'f_Bmag', 'e_Vmag', 'f_Vmag', 'gmag', 'e_gmag', 'f_gmag', 'rmag', 'e_rmag', 'f_rmag', 'imag', 'e_imag', 'f_imag', 'LEDA', 'Hmag'])

	Vizier.ROW_LIMIT = -1
	v.ROW_LIMIT = -1

	# note: height and width seem to be mixed up in some astropy versions
	result = v.query_region(coord.ICRS(ra=st, dec=sde, unit=(u.deg, u.deg)), width=str(2*rad)+"d", height=str(step)+"d", catalog='I/322A')
	#width=str(2*rad)+"d", height=str(step)+"d", catalog='I/322A')

	table=result[0]
#	print(table.colnames)
#	pdb.set_trace()	

	if ix == 0:
		ra = np.array(table[:]['_RAJ2000'])
		dec = np.array(table[:]['_DEJ2000'])
		ucac = np.array(table[:]['UCAC4'])
		pmra = np.array(table[:]['pmRA'])
		e_pmra = np.array(table[:]['e_pmRA'])
		pmde = np.array(table[:]['pmDE'])
		e_pmde =  np.array(table[:]['e_pmDE'])
		tycho = np.array(table[:]['Tycho-2'])
		mkey =  np.array(table[:]['_2Mkey'])
		jmag = np.array(table[:]['Jmag'])
		e_jmag = np.array(table[:]['e_Jmag'])
		hmag = np.array(table[:]['Hmag'])
		e_hmag = np.array(table[:]['e_Hmag'])
		kmag = np.array(table[:]['Kmag'])
		e_kmag = np.array(table[:]['e_Kmag'])
		bmag = np.array(table[:]['Bmag'])
		e_bmag = np.array(table[:]['e_Bmag'])
		f_bmag = np.array(table[:]['f_Bmag'])
		vmag =  np.array(table[:]['Vmag'])
		e_vmag = np.array(table[:]['e_Vmag'])
		f_vmag = np.array(table[:]['f_Vmag'])
		gmag = np.array(table[:]['gmag'])
		e_gmag = np.array(table[:]['e_gmag'])
		f_gmag = np.array(table[:]['f_gmag'])
		rmag = np.array(table[:]['rmag'])
		e_rmag = np.array(table[:]['e_rmag'])
		f_rmag =np.array(table[:]['f_rmag'])
		imag = np.array(table[:]['imag'])
		e_imag = np.array(table[:]['e_imag'])
		f_imag = np.array(table[:]['f_imag'])
		leda = np.array(table[:]['LEDA'])

	else:
		ra = np.concatenate((ra,np.array(table[:]['_RAJ2000'])))
		dec = np.concatenate((dec,np.array(table[:]['_DEJ2000'])))
		ucac = np.concatenate((ucac,np.array(table[:]['UCAC4'])))
		pmra = np.concatenate((pmra,np.array(table[:]['pmRA'])))
		e_pmra = np.concatenate((e_pmra,np.array(table[:]['e_pmRA'])))
		pmde = np.concatenate((pmde,np.array(table[:]['pmDE'])))
		e_pmde = np.concatenate((e_pmde,np.array(table[:]['e_pmDE'])))
		tycho = np.concatenate((tycho,np.array(table[:]['Tycho-2'])))
		mkey = np.concatenate((mkey,np.array(table[:]['_2Mkey'])))
		jmag = np.concatenate((jmag,np.array(table[:]['Jmag'])))
		e_jmag = np.concatenate((e_jmag,np.array(table[:]['e_Jmag'])))
		hmag = np.concatenate((hmag,np.array(table[:]['Hmag'])))
		e_hmag = np.concatenate((e_hmag,np.array(table[:]['e_Hmag'])))
		kmag = np.concatenate((kmag,np.array(table[:]['Kmag'])))
		e_kmag = np.concatenate((e_kmag,np.array(table[:]['e_Kmag'])))
		bmag = np.concatenate((bmag,np.array(table[:]['Bmag'])))
		e_bmag = np.concatenate((e_bmag,np.array(table[:]['e_Bmag'])))
		f_bmag = np.concatenate((f_bmag,np.array(table[:]['f_Bmag'])))
		vmag = np.concatenate((vmag,np.array(table[:]['Vmag'])))
		e_vmag = np.concatenate((e_vmag,np.array(table[:]['e_Vmag'])))
		f_vmag = np.concatenate((f_vmag,np.array(table[:]['f_Vmag'])))
		gmag = np.concatenate((gmag,np.array(table[:]['gmag'])))
		e_gmag = np.concatenate((e_gmag,np.array(table[:]['e_gmag'])))
		f_gmag = np.concatenate((f_gmag,np.array(table[:]['f_gmag'])))	
		rmag = np.concatenate((rmag,np.array(table[:]['rmag'])))
		e_rmag = np.concatenate((e_rmag,np.array(table[:]['e_rmag'])))
		f_rmag = np.concatenate((f_rmag,np.array(table[:]['f_rmag'])))		
		imag = np.concatenate((imag,np.array(table[:]['imag'])))
		e_imag = np.concatenate((e_imag,np.array(table[:]['e_imag'])))
		f_imag = np.concatenate((f_imag,np.array(table[:]['f_imag'])))
		leda = np.concatenate((leda,np.array(table[:]['LEDA'])))
	
	st = st + step
	ix = ix+1
			

# throw out duplicates, and restrict to search radius
unq,unq_ix=np.unique(ucac,return_index=True)	
dis = np.sqrt( np.power((ra[unq_ix]-sra),2) + np.power((dec[unq_ix]-sde),2) )
keep = np.where(dis < rad)

#pylab.ion()
#plt.plot(ra[unq_ix[keep[0]]],dec[unq_ix[keep[0]]],'b.')
#pylab.draw()

print('Writing UCAC fits file')
cols = fits.ColDefs([\
fits.Column(name='_RAJ2000', format='D', array=ra[unq_ix[keep[0]]]),\
fits.Column(name='_DEJ2000', format='D', array=dec[unq_ix[keep[0]]]),\
fits.Column(name='UCAC4', format='20A', array=ucac[unq_ix[keep[0]]]),\
fits.Column(name='PMRA', format='D', array=pmra[unq_ix[keep[0]]]),\
fits.Column(name='E_PMRA', format='D', array=e_pmra[unq_ix[keep[0]]]),\
fits.Column(name='PMDE', format='D', array=pmde[unq_ix[keep[0]]]),\
fits.Column(name='E_PMDE', format='D', array=e_pmde[unq_ix[keep[0]]]),\
fits.Column(name='TYCHO_2', format='20A', array=tycho[unq_ix[keep[0]]]),\
fits.Column(name='_2MKEY', format='D', array=mkey[unq_ix[keep[0]]]),\
fits.Column(name='JMAG', format='D', array=jmag[unq_ix[keep[0]]]),\
fits.Column(name='E_JMAG', format='D', array=e_jmag[unq_ix[keep[0]]]),\
fits.Column(name='HMAG', format='D', array=hmag[unq_ix[keep[0]]]),\
fits.Column(name='E_HMAG', format='D', array=e_hmag[unq_ix[keep[0]]]),\
fits.Column(name='KMAG', format='D', array=kmag[unq_ix[keep[0]]]),\
fits.Column(name='E_KMAG', format='D', array=e_kmag[unq_ix[keep[0]]]),\
fits.Column(name='BMAG', format='D', array=bmag[unq_ix[keep[0]]]),\
fits.Column(name='E_BMAG', format='D', array=e_bmag[unq_ix[keep[0]]]),\
fits.Column(name='F_BMAG', format='A', array=f_bmag[unq_ix[keep[0]]]),\
fits.Column(name='VMAG', format='D', array=vmag[unq_ix[keep[0]]]),\
fits.Column(name='E_VMAG', format='D', array=e_vmag[unq_ix[keep[0]]]),\
fits.Column(name='F_VMAG', format='A', array=f_vmag[unq_ix[keep[0]]]),\
fits.Column(name='GMAG', format='D', array=gmag[unq_ix[keep[0]]]),\
fits.Column(name='E_GMAG', format='D', array=e_gmag[unq_ix[keep[0]]]),\
fits.Column(name='F_GMAG', format='A', array=f_gmag[unq_ix[keep[0]]]),\
fits.Column(name='RMAG', format='D', array=rmag[unq_ix[keep[0]]]),\
fits.Column(name='E_RMAG', format='D', array=e_rmag[unq_ix[keep[0]]]),\
fits.Column(name='F_RMAG', format='A', array=f_rmag[unq_ix[keep[0]]]),\
fits.Column(name='IMAG', format='D', array=imag[unq_ix[keep[0]]]),\
fits.Column(name='E_IMAG', format='D', array=e_imag[unq_ix[keep[0]]]),\
fits.Column(name='F_IMAG', format='A', array=f_imag[unq_ix[keep[0]]]),\
fits.Column(name='LEDA', format='D', array=leda[unq_ix[keep[0]]])])

tbhdu = fits.new_table(cols)
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('catalogs/ucac.fits',clobber=True)

#pdb.set_trace()	



##########################################################
###### 2MASS
	
st = (sra-rad)
ix=0

while st <= sra+rad:

	print('Querying 2MASS, tile '+str(ix+1)+' of '+str(ntiles))
	t0 = time.time()
	#v = Vizier(columns=['_RAJ2000', '_DEJ2000','2MASS', 'Jmag', 'e_Jmag',\
	#'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag', 'Qflg', 'Rflg', 'Bflg', 'Cflg', 'Xflg', 'Aflg', \
	#'prox','Cntr'])
	
        #v = Vizier(columns=[':'])
	# only query for those columns missing from the default output, which are returned anyway
	# same is implemented for NOMAD+SDSS
     #   v = Vizier(columns=['prox','Cntr'])
	v = Vizier(columns=['all'])

	Vizier.ROW_LIMIT = -1
	v.ROW_LIMIT = -1
	result = v.query_region(coord.ICRS(ra=st, dec=sde, unit=(u.deg, u.deg)), width=str(2*rad)+"d", height=str(step)+"d", catalog='II/246')

	#result = v.query_region(coord.ICRS(ra=sra, dec=sde, unit=(u.deg, u.deg)), \
	#radius=Angle(rad, "deg"), catalog='II/246')

	#width=str(2*rad)+"d", height=str(step)+"d", catalog='II/246')

	#str(2*rad)+

	table=result[0]
	#print(table.colnames)
	#pbd.set_trace()
	#x=np.array(table[:]['Plx'])
	#y=np.array(table[:]['e_Plx'])

	if ix == 0:
		ra = np.array(table[:]['_RAJ2000'])
		dec = np.array(table[:]['_DEJ2000'])
		mass = np.array(table[:]['_2MASS'])
		jmag = np.array(table[:]['Jmag'])
		e_jmag = np.array(table[:]['e_Jmag'])
		hmag = np.array(table[:]['Hmag'])
		e_hmag = np.array(table[:]['e_Hmag'])
		kmag = np.array(table[:]['Kmag'])
		e_kmag = np.array(table[:]['e_Kmag'])
		qflg = np.array(table[:]['Qflg'])
		rflg = np.array(table[:]['Rflg'])
		bflg = np.array(table[:]['Bflg'])
		cflg = np.array(table[:]['Cflg'])
		xflg = np.array(table[:]['Xflg'])
		aflg = np.array(table[:]['Aflg'])
		prox = np.array(table[:]['prox'])
		centr = np.array(table[:]['Cntr'])
				
	else:
		ra = np.concatenate((ra,np.array(table[:]['_RAJ2000'])))
		dec = np.concatenate((dec,np.array(table[:]['_DEJ2000'])))
		mass = np.concatenate((mass,np.array(table[:]['_2MASS'])))
		jmag = np.concatenate((jmag,np.array(table[:]['Jmag'])))
		e_jmag = np.concatenate((e_jmag,np.array(table[:]['e_Jmag'])))
		hmag = np.concatenate((hmag,np.array(table[:]['Hmag'])))
		e_hmag = np.concatenate((e_hmag,np.array(table[:]['e_Hmag'])))
		kmag = np.concatenate((kmag,np.array(table[:]['Kmag'])))
		e_kmag = np.concatenate((e_kmag,np.array(table[:]['e_Kmag'])))
		qflg = np.concatenate((qflg,np.array(table[:]['Qflg'])))
		rflg = np.concatenate((rflg,np.array(table[:]['Rflg'])))
		bflg = np.concatenate((bflg,np.array(table[:]['Bflg'])))
		cflg = np.concatenate((cflg,np.array(table[:]['Cflg'])))
		xflg = np.concatenate((xflg,np.array(table[:]['Xflg'])))
		aflg = np.concatenate((aflg,np.array(table[:]['Aflg'])))
		prox = np.concatenate((prox,np.array(table[:]['prox'])))		
		centr = np.concatenate((centr,np.array(table[:]['Cntr'])))
		
	st = st + step
	ix = ix+1

# throw out duplicates, and restrict to search radius
unq,unq_ix=np.unique(centr,return_index=True)	
dis = np.sqrt( np.power((ra[unq_ix]-sra),2) + np.power((dec[unq_ix]-sde),2) )
keep = np.where(dis < rad)

#pylab.ion()
#plt.plot(ra[unq_ix[keep[0]]],dec[unq_ix[keep[0]]],'b.')
#pylab.draw()

print('Writing 2MASS fits file')
cols = fits.ColDefs([\
fits.Column(name='_RAJ2000', format='D', array=ra[unq_ix[keep[0]]]),\
fits.Column(name='_DEJ2000', format='D', array=dec[unq_ix[keep[0]]]),\
fits.Column(name='_2MASS', format='30A', array=mass[unq_ix[keep[0]]]),\
fits.Column(name='JMAG', format='D', array=jmag[unq_ix[keep[0]]]),\
fits.Column(name='E_JMAG', format='D', array=e_jmag[unq_ix[keep[0]]]),\
fits.Column(name='HMAG', format='D', array=hmag[unq_ix[keep[0]]]),\
fits.Column(name='E_HMAG', format='D', array=e_hmag[unq_ix[keep[0]]]),\
fits.Column(name='KMAG', format='D', array=kmag[unq_ix[keep[0]]]),\
fits.Column(name='E_KMAG', format='D', array=e_kmag[unq_ix[keep[0]]]),\
fits.Column(name='QFLG', format='3A', array=qflg[unq_ix[keep[0]]]),\
fits.Column(name='RFLG', format='3A', array=rflg[unq_ix[keep[0]]]),\
fits.Column(name='BFLG', format='3A', array=bflg[unq_ix[keep[0]]]),\
fits.Column(name='CFLG', format='3A', array=cflg[unq_ix[keep[0]]]),\
fits.Column(name='XFLG', format='3A', array=xflg[unq_ix[keep[0]]]),\
fits.Column(name='AFLG', format='3A', array=aflg[unq_ix[keep[0]]]),\
fits.Column(name='prox', format='D', array=prox[unq_ix[keep[0]]]),\
fits.Column(name='CNTR', format='D', array=centr[unq_ix[keep[0]]])])

tbhdu = fits.new_table(cols)
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('catalogs/2mass.fits',clobber=True)
#print(time.time()-t0)

#pdb.set_trace()


##########################################################
###### NOMAD
	
st = (sra-rad)
ix=0

while st <= sra+rad:

	print('Querying NOMAD, tile '+str(ix+1)+' of '+str(ntiles))
	t0 = time.time()
	
	
	#v = Vizier(columns=['YM'])
	v = Vizier(columns=['all'])
	Vizier.ROW_LIMIT = -1
	v.ROW_LIMIT = -1
	result = v.query_region(coord.ICRS(ra=st, dec=sde, unit=(u.deg, u.deg)), width=str(2*rad)+"d", height=str(step)+"d", catalog='I/297')

	#width=str(2*rad)+"d", height=str(step)+"d", catalog='I/297')

	#str(2*rad)+

	table=result[0]
	#print(table.colnames)
	#pbd.set_trace()
	#x=np.array(table[:]['Plx'])
	#y=np.array(table[:]['e_Plx'])

	if ix == 0:
		ra = np.array(table[:]['_RAJ2000'])
		dec = np.array(table[:]['_DEJ2000'])
		nomad = np.array(table[:]['NOMAD1'])
		jmag = np.array(table[:]['Jmag'])
		hmag = np.array(table[:]['Hmag'])
		kmag = np.array(table[:]['Kmag'])
		pmra = np.array(table[:]['pmRA'])
		e_pmra = np.array(table[:]['e_pmRA'])
		pmde = np.array(table[:]['pmDE'])
		e_pmde =  np.array(table[:]['e_pmDE'])
		ym =  np.array(table[:]['YM'])
		
				
	else:
		ra = np.concatenate((ra,np.array(table[:]['_RAJ2000'])))
		dec = np.concatenate((dec,np.array(table[:]['_DEJ2000'])))
		nomad = np.concatenate((nomad,np.array(table[:]['NOMAD1'])))
		jmag = np.concatenate((jmag,np.array(table[:]['Jmag'])))
		hmag = np.concatenate((hmag,np.array(table[:]['Hmag'])))
		kmag = np.concatenate((kmag,np.array(table[:]['Kmag'])))
		pmra = np.concatenate((pmra,np.array(table[:]['pmRA'])))
		e_pmra = np.concatenate((e_pmra,np.array(table[:]['e_pmRA'])))
		pmde = np.concatenate((pmde,np.array(table[:]['pmDE'])))
		e_pmde = np.concatenate((e_pmde,np.array(table[:]['e_pmDE'])))
		ym = np.concatenate((ym,np.array(table[:]['YM'])))

		
	st = st + step
	ix = ix+1

# throw out duplicates, and restrict to search radius
unq,unq_ix=np.unique(nomad,return_index=True)	
dis = np.sqrt( np.power((ra[unq_ix]-sra),2) + np.power((dec[unq_ix]-sde),2) )
keep = np.where(dis < rad)

#pylab.ion()
#plt.plot(ra[unq_ix[keep[0]]],dec[unq_ix[keep[0]]],'b.')
#pylab.draw()

print('Writing NOMAD fits file')
cols = fits.ColDefs([\
fits.Column(name='_RAJ2000', format='D', array=ra[unq_ix[keep[0]]]),\
fits.Column(name='_DEJ2000', format='D', array=dec[unq_ix[keep[0]]]),\
fits.Column(name='NOMAD1', format='30A', array=nomad[unq_ix[keep[0]]]),\
fits.Column(name='JMAG', format='D', array=jmag[unq_ix[keep[0]]]),\
fits.Column(name='HMAG', format='D', array=hmag[unq_ix[keep[0]]]),\
fits.Column(name='KMAG', format='D', array=kmag[unq_ix[keep[0]]]),\
fits.Column(name='YM', format='5A', array=ym[unq_ix[keep[0]]]),\
fits.Column(name='PMRA', format='D', array=pmra[unq_ix[keep[0]]]),\
fits.Column(name='E_PMRA', format='D', array=e_pmra[unq_ix[keep[0]]]),\
fits.Column(name='PMDE', format='D', array=pmde[unq_ix[keep[0]]]),\
fits.Column(name='E_PMDE', format='D', array=e_pmde[unq_ix[keep[0]]])])

tbhdu = fits.new_table(cols)
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('catalogs/nomad.fits',clobber=True)
#print(time.time()-t0)

#pdb.set_trace()


##########################################################
###### SDSS

st = (sra-rad)
ix=0

while st <= sra+rad:
	
	print('Querying SDSS, tile '+str(ix+1)+' of '+str(ntiles))

#    v = Vizier(columns=['SDSS-ID', 'objID', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE'], column_filters={"rmag":"<20"})
#	v = Vizier(columns=['SDSS-ID', 'objID', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE'],)
	v = Vizier(columns=['all'], column_filters={"rmag":"<20"})
	Vizier.ROW_LIMIT = -1
	v.ROW_LIMIT = -1
	result = v.query_region(coord.ICRS(ra=st, dec=sde, unit=(u.deg, u.deg)), width=str(2*rad)+"d", height=str(step)+"d", catalog='V/139')
		
	# Sloan isn't all-sky, so skip if there are no stars queried
	if len(result) == 0:
		st = st + step
		ix = ix+1
		continue
					
	table=result[0]
	#print(table.colnames)
	
	if ix == 0:
		ra = np.array(table[:]['_RAJ2000'])
		dec = np.array(table[:]['_DEJ2000'])
		qmode = np.array(table[:]['q_mode'])
		cl =  np.array(table[:]['cl'])
		sdss9 =  np.array(table[:]['SDSS9'])
		sdssid =  np.array(table[:]['SDSS-ID'])
		objid =  np.array(table[:]['objID'])
		q =  np.array(table[:]['Q'])
		umag = np.array(table[:]['umag'])
		e_umag = np.array(table[:]['e_umag'])
		gmag = np.array(table[:]['gmag'])
		e_gmag = np.array(table[:]['e_gmag'])
		rmag = np.array(table[:]['rmag'])
		e_rmag = np.array(table[:]['e_rmag'])
		imag = np.array(table[:]['imag'])
		e_imag = np.array(table[:]['e_imag'])
		zmag = np.array(table[:]['zmag'])
		e_zmag = np.array(table[:]['e_zmag'])		
		pmra = np.array(table[:]['pmRA'])
		e_pmra = np.array(table[:]['e_pmRA'])
		pmde = np.array(table[:]['pmDE'])
		e_pmde =  np.array(table[:]['e_pmDE'])
		
	else:
		ra = np.concatenate((ra,np.array(table[:]['_RAJ2000'])))
		dec = np.concatenate((dec,np.array(table[:]['_DEJ2000'])))
		qmode = np.concatenate((qmode,np.array(table[:]['q_mode'])))
		cl = np.concatenate((cl,np.array(table[:]['cl'])))
		sdss9 = np.concatenate((sdss9,np.array(table[:]['SDSS9'])))
		sdssid = np.concatenate((sdssid,np.array(table[:]['SDSS-ID'])))
		objid = np.concatenate((objid,np.array(table[:]['objID'])))
		q = np.concatenate((q,np.array(table[:]['Q'])))
		umag = np.concatenate((umag,np.array(table[:]['umag'])))
		e_umag = np.concatenate((e_umag,np.array(table[:]['e_umag'])))
		gmag = np.concatenate((gmag,np.array(table[:]['gmag'])))
		e_gmag = np.concatenate((e_gmag,np.array(table[:]['e_gmag'])))
		rmag = np.concatenate((rmag,np.array(table[:]['rmag'])))
		e_rmag = np.concatenate((e_rmag,np.array(table[:]['e_rmag'])))
		imag = np.concatenate((imag,np.array(table[:]['imag'])))
		e_imag = np.concatenate((e_imag,np.array(table[:]['e_imag'])))
		zmag = np.concatenate((zmag,np.array(table[:]['zmag'])))
		e_zmag = np.concatenate((e_zmag,np.array(table[:]['e_zmag'])))
		pmra = np.concatenate((pmra,np.array(table[:]['pmRA'])))
		e_pmra = np.concatenate((e_pmra,np.array(table[:]['e_pmRA'])))
		pmde = np.concatenate((pmde,np.array(table[:]['pmDE'])))
		e_pmde = np.concatenate((e_pmde,np.array(table[:]['e_pmDE'])))

	st = st + step
	ix = ix+1


# throw out duplicates, and restrict to search radius
unq,unq_ix=np.unique(sdss9,return_index=True)	
dis = np.sqrt( np.power((ra[unq_ix]-sra),2) + np.power((dec[unq_ix]-sde),2) )
keep = np.where(dis < rad)

#pylab.ion()
#plt.plot(ra[unq_ix[keep[0]]],dec[unq_ix[keep[0]]],'b.')
#pylab.draw()

print('Writing SDSS fits file')
cols = fits.ColDefs([\
fits.Column(name='_RAJ2000', format='D', array=ra[unq_ix[keep[0]]]),\
fits.Column(name='_DEJ2000', format='D', array=dec[unq_ix[keep[0]]]),\
fits.Column(name='Q_MODE', format='A', array=qmode[unq_ix[keep[0]]]),\
fits.Column(name='CL', format='I', array=cl[unq_ix[keep[0]]]),\
fits.Column(name='SDSS9', format='30A', array=sdss9[unq_ix[keep[0]]]),\
fits.Column(name='SDSS_ID', format='30A', array=sdssid[unq_ix[keep[0]]]),\
fits.Column(name='OBJID', format='I', array=objid[unq_ix[keep[0]]]),\
fits.Column(name='Q', format='A', array=q[unq_ix[keep[0]]]),\
fits.Column(name='UMAG', format='D', array=umag[unq_ix[keep[0]]]),\
fits.Column(name='E_UMAG', format='D', array=e_umag[unq_ix[keep[0]]]),\
fits.Column(name='GMAG', format='D', array=gmag[unq_ix[keep[0]]]),\
fits.Column(name='E_GMAG', format='D', array=e_gmag[unq_ix[keep[0]]]),\
fits.Column(name='RMAG', format='D', array=rmag[unq_ix[keep[0]]]),\
fits.Column(name='E_RMAG', format='D', array=e_rmag[unq_ix[keep[0]]]),\
fits.Column(name='IMAG', format='D', array=imag[unq_ix[keep[0]]]),\
fits.Column(name='E_IMAG', format='D', array=e_imag[unq_ix[keep[0]]]),\
fits.Column(name='ZMAG', format='D', array=zmag[unq_ix[keep[0]]]),\
fits.Column(name='E_ZMAG', format='D', array=e_zmag[unq_ix[keep[0]]]),\
fits.Column(name='PMRA', format='D', array=pmra[unq_ix[keep[0]]]),\
fits.Column(name='E_PMRA', format='D', array=e_pmra[unq_ix[keep[0]]]),\
fits.Column(name='PMDE', format='D', array=pmde[unq_ix[keep[0]]]),\
fits.Column(name='E_PMDE', format='D', array=e_pmde[unq_ix[keep[0]]])])

tbhdu = fits.new_table(cols)
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('catalogs/sdss.fits',clobber=True)
#print(time.time()-t0)

#pdb.set_trace()	







	
	
	

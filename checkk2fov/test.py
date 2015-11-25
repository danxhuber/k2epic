import pdb
import numpy as np
import matplotlib.pyplot as plt


f = open('/Users/daniel/science/K2/EPIC/deliveries/d14260_03_epic_c7_dmc.mrg', 'r')
n=0
for line in f:
    n=n+1
    print(n)
    
f.close()
ra=np.zeros(n)
dec=np.zeros(n)

f = open('/Users/daniel/science/K2/EPIC/deliveries/d14260_03_epic_c7_dmc.mrg', 'r')
n=0
for line in f:

	if (n % 10000 != 0):
		n=n+1
		continue

	print(n)
	line = line.strip()
	columns = line.split("|")
	ra[n]=columns[9]
	dec[n]=columns[10]
	n=n+1
	
use=np.where(ra > 0)[0]
ra=ra[use]
dec=dec[use]




    
pdb.set_trace()

import sys
import numpy as np
import re

bf = sys.argv[1]
nbins = int(sys.argv[2])
field = int(sys.argv[3])


fn = bf 
f = open(fn,'r')
datar = f.read().split('\n') 

data = []
for j in range(len(datar)):
	lsp = re.split(' |\t',datar[j])
	try:
		data.insert(j,float(lsp[field]))
	except:
		pass

thist,bedge = np.histogram(data,nbins)

l = len(thist)
for j in range(l):
	print '%6.3e %6.3e' % (bedge[j],thist[j])



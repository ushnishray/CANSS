
import sys
import numpy as np

bf = sys.argv[1]
nbins = int(sys.argv[2])

fn = bf 
f = open(fn,'r')
datar = f.read().split('\n') 

data = []
for j in range(len(datar)):
	try:
		data.insert(j,float(datar[j]))
	except:
		pass

thist,bedge = np.histogram(data,nbins)

l = len(thist)
for j in range(l):
	print '%6.3e %6.3e' % (bedge[j],thist[j])




import sys
import numpy as np
import re

if(len(sys.argv)<2):
	print "Usage is: <no. of files> <base file name> <no. of bins> [<field no.>]"
	sys.exit()

nf = int(sys.argv[1])
bf = sys.argv[2]
nbins = int(sys.argv[3])
if(len(sys.argv)<4):
	field = 0
else:
	field = int(sys.argv[4])

bedge = []
thist = []
for i in range(0,nf):

	fn = bf + str(i)
	f = open(fn,'r')
	datar = f.read().split('\n') 
	
	data = []
	for j in range(len(datar)):
		lsp = re.split(' |\t',datar[j])
		try:
			data.insert(j,float(lsp[field]))
		except:
			print lsp
			pass

	if(i == 0):
		thist,bedge = np.histogram(data,nbins)
	else:
		hist,ledge = np.histogram(data,bedge)
		thist = [x + y for x,y in zip(thist,hist)]

fhist = [float(x)/nf for x in thist]

total = np.sum(fhist)
l = len(fhist)
for j in range(l):
	print '%6.3e %6.3e' % (bedge[j],fhist[j]/total)



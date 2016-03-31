import numpy
import sys
 
def acf(series):
	n = len(series)
	data = numpy.asarray(series)
	mean = numpy.mean(data)
	c0 = numpy.sum((data-mean)**2)/float(n)

	def  r(h):
		acf_lag = ((data[:n-h] - mean)*(data[h:] - mean)).sum() / float(n) / c0
		return round(acf_lag,3)

	x = numpy.arange(n)
	acf_coeffs = map(r,x)
	return acf_coeffs

datf = sys.argv[1]
field = int(sys.argv[2])
ns = int(sys.argv[3])

f = open(datf,'r')
x = []
x2 = []

binc = 0
n = 0
l = []
for a in f:
	d = float(a.split(' ')[field])
	l.append(d)
	n += 1
	if(n%ns == 0):
		af = acf(l)
		p = 0
		for j in af:
			if(binc == 0):
				x.append(j)
				x2.append(j*j)
			else:
				x[p] += j
				x2[p] += j*j
			p += 1	
		n = 0	
		l = []
		binc += 1

f.close()

for i in range(ns):
	x[i] /= binc 
	x2[i] /= binc
	x2[i] = ((x2[i]-x[i]*x[i])/binc)**0.5

	print i,x[i],x2[i]


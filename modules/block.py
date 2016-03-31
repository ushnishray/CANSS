import numpy as np
import sys

def block(data, start, piece, tsize, level, x, x2):

	if(piece == 1):
		try:
			x[level] += data[start]
			x2[level] += data[start]**2
#			print level,data[start],x[level]
		except:
			print level

		return data[start]

	a = block(data,start,piece/2,tsize, level+1, x, x2)
	b = block(data,start+piece/2,piece/2,tsize, level+1, x, x2)
	val = (a+b)/2
	x[level] += val
	x2[level] += val**2
	
#	print level,val,x[level]
	return val

def test():
	data = [4.,5.,6.,7.]
	x = []
	x2 = []
	for i in range(3):
		x.append(0)
		x2.append(0)

	a = block(data,0,4,4,2,x,x2)

	tsize = 4
	r = 3
	for i in range(r-1):
		ns = tsize/2**i
		x[i] /= tsize
		x2[i] /= ns
		error = ((x2[i]-x[i]**2)/(ns-1))**0.5
		unc = error/(2*(ns-1))**0.5	
		print i,x[i],error,unc

def regular():
	df = open(sys.argv[1],'r')
	field = int(sys.argv[2])
	data = []

	for dl in df:
		d = dl.split(' ')
		data.append(float(d[field]))
	df.close()

	l = len(data)
	r = int(np.floor(np.log(l)/np.log(2)))

	data = data[-2**r:]
	x = []
	x2 = []
	for i in range(r+1):
		x.append(0)
		x2.append(0)

	tsize = 2**r
	a = block(data,0,tsize,tsize,0,x,x2)

	print "Readjusted data size: ",len(data)," Levels: ",r
#	for i in range(r+1):
#		print x[i],x2[i]

	for i in range(1,r+1):

		ns = 2**i
		x[i] /= ns
 		x2[i] /= ns
		error = ((x2[i]-x[i]*x[i])/(ns-1))**0.5
		unc = error/(2*(ns-1))**0.5	
		print i,x[i],error,unc

regular()

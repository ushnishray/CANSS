import os
import sys
import math

if(len(sys.argv)<8):
	print 'format is: <beta file> <base hist file> <no. of hist files> <width> <start val> <end val> <no. of walkers>\n'
	sys.exit()

fin = open(sys.argv[1])
vbetas = fin.read().split('\n')
fin.close()
betas = []
for beta in vbetas:
	try:
		float(beta)
		betas.append(beta)
	except:
		pass


baseHist = sys.argv[2]
nhist = int(sys.argv[3])

width = float(sys.argv[4])
sb = float(sys.argv[5])
eb = float(sys.argv[6])

nw = int(sys.argv[7])
norm = 1.0/(nw*len(betas));

didx = []
dpairs = []
x = sb
li = 0
while x<eb:
	didx.append(x)
	dpairs.append((li,0))

	li += 1
	x += width
didx.append(x)
dpairs.append((li,0))

dgather = dict(dpairs)
dgather2 = dict(dpairs)

for r in range(0,nhist):

	dgatherL = dict(dpairs)

	i = 0
	for beta in betas:
		float(beta)
		sdir = 'set.' + str(i) + '/'
		ffile = baseHist + str(r)
		fin = open(sdir + ffile)
		for qvall in fin:
			qval = float(qvall.split()[0]) - sb

			bdict = math.floor((qval)/width)
			try:
				dgatherL[bdict] += 1
			except:
				dgatherL[li] += 1
		fin.close()


	#normalize
	lnorm = 0.0
	for key in dgather.iterkeys():
		lnorm += dgatherL[key]
	
	lnorm = 1.0/lnorm

	print dgatherL

	#We have the bins
	#print dgatherL
	for key in dgather.iterkeys():
		try:
			val = dgatherL[key]*lnorm
			dgather[key] += val
			dgather2[key] += val*val
		except:
			pass
	
		
for key in dgather.iterkeys():
	dgather[key] /= nhist
	dgather2[key] /= nhist
	dgather2[key] = ((dgather2[key]-dgather[key]**2)/(nhist-1))**0.5 
	
	print '% 10.6e\t% 10.6e\t% 10.6e' % (didx[key], dgather[key], dgather2[key])



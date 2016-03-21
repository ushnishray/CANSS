import os
import sys

fin = open(sys.argv[1])
datf = sys.argv[2]

betas = fin.read().split('\n')
fin.close()

w = ''
i = 0
for beta in betas:
	try:
		lb = float(beta)
		sdir = 'set.' + str(i)


		psi = 0.0
		psi2 = 0.0
		q = 0.0
		q2 = 0.0
		v = 0.0
		v2 = 0.0
		h = h2 = 0.0
		n = 0

		dfo = open(sdir + '/' + datf,'r')	
		for l in dfo:		
			st = l.split()			
			psi += float(st[1])
			psi2 += float(st[1])**2
			q += float(st[2])
			q2 += float(st[2])**2
			v += float(st[5])
			v2 += float(st[5])**2
			
			lh = lb*float(st[2]) - float(st[1])	
			h += lh
			h2 += lh**2
			
			n += 1	
		dfo.close()

		psi /= n; psi2 /= n
		psi2 = ((psi2 - psi**2)/(n-1))**0.5
		q /= n; q2 /= n
		q2 = ((q2-q**2)/(n-1))**0.5
		v /=n; v2 /= n
		v2 = ((v2-v**2)/(n-1))**0.5
		h /=n; h2 /= n
		h2 = ((h2-h**2)/(n-1))**0.5

		w += beta + '\t' + str(psi) + '\t' + str(psi2) + '\t' + str(q) + '\t' + str(q2) + '\t' + str(v) + '\t' + str(v2) + '\t' + str(h) + '\t' + str(h2) + '\n'
		i += 1
	except:
		pass		


print w

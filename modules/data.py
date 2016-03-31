import os
import sys

#offset = 0
offset = 15

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
		psie = 0.0
		psie2 = 0.0

		q = 0.0
		q2 = 0.0
		qe = 0.0
		qe2 = 0.0

		v = 0.0
		v2 = 0.0
		ve = 0.0
		ve2 = 0.0

		h = h2 = 0.0
		n = 0

		dfo = open(sdir + '/' + datf,'r')	
		for l in dfo:		
			st = l.split()			
			psi += float(st[offset+1])
			psi2 += float(st[offset+1])**2
			psie += float(st[offset+2])
			psie2 += float(st[offset+2])**2
			
			q += float(st[3+offset])
			q2 += float(st[3+offset])**2
			qe += float(st[4+offset])
			qe2 += float(st[4+offset])**2

			v += float(st[9+offset])
			v2 += float(st[9+offset])**2			
			ve += float(st[10+offset])
			ve2 += float(st[10+offset])**2			

			lh = lb*float(st[3+offset]) - float(st[1+offset])	
			h += lh
			h2 += lh**2
			
			n += 1	
		dfo.close()

		psi /= n; psi2 /= n
		psi2 = ((psi2 - psi**2)/(n))**0.5
		psie /= n; psie2 /= n
		psie2 = ((psie2 - psie**2)/(n))**0.5

		q /= n; q2 /= n
		q2 = ((q2-q**2)/(n))**0.5
		qe /= n; qe2 /= n
		qe2 = ((qe2-qe**2)/(n))**0.5

		v /=n; v2 /= n
		v2 = ((v2-v**2)/(n))**0.5
		ve /=n; ve2 /= n
		ve2 = ((ve2-ve**2)/(n))**0.5

		h /=n; h2 /= n
		h2 = ((h2-h**2)/(n))**0.5

#		w += beta + '\t' + str(psi) + '\t' + str(psi2) + '\t' + str(q) + '\t' + str(q2) + '\t' + str(v) + '\t' + str(v2) + '\t' + str(h) + '\t' + str(h2) + '\n'
		w += beta + '\t' + str(psi) + '\t' + str(psie) + '\t' + str(q) + '\t' + str(qe) + '\t' + str(v) + '\t' + str(ve) + '\t' + str(h) + '\t' + str(h2) + '\t'
		w += str(psi2) + '\t' + str(psie2) + '\t' + str(q2) + '\t' + str(qe2) + '\t' + str(v2) + '\t' + str(ve2) + '\n'
		i += 1
	except:
		pass		


print w

import os
import sys

fin = open(sys.argv[1])
betas = fin.read().split('\n')
fin.close()

w = ''
bexec = '/home/ushnish/research/Development/checkouts/DMCsns/cppmods/models/DBM/d1/DEX.branch'
i = 0

nproc = 2
tproc = 8
for beta in betas:
	try:
		float(beta)
		sdir = 'set.' + str(i)
		os.system('mkdir ' + sdir)
		os.system('cp -r params ' + sdir + '/')

		w += 'cd ../' + sdir + '\n'
		w += 'mpirun -np ' + str(nproc) + ' ' + bexec + ' ' + './params/param.txt ' + str(beta) + ' >> out &\n\n'
		i += 1
		if((i*nproc)%tproc == 0):
			w += 'wait\n'
	except:
		pass		
w += 'wait\n'
print w

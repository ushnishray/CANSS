import os
import sys

fin = open(sys.argv[1])
betas = fin.read().split('\n')
fin.close()

w = ''
bexec = '/home/ushnish/research/dmcsns/models/DBM/d1/DEX.nb' 
#bexec = '/home/ushnish/research/dmcsns/models/SNS/d1b/DEX.branch'
i = 0

nproc = 31
tproc = 32
for beta in betas:
	try:
		float(beta)
		sdir = 'set.' + str(i)
		os.system('mkdir ' + sdir)
		os.system('cp -r params ' + sdir + '/')

		w += 'cd ../' + sdir + '\n'
		w += 'mpirun -np ' + str(nproc) + ' ' + bexec + ' ' + './params/param.txt ' + str(beta) + ' >> out \n\n'
		i += 1
		if((i*nproc)%tproc == 0):
			w += 'wait\n'
	except:
		pass		
w += 'wait\n'
print w

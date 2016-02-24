import os
import sys

fin = open(sys.argv[1])
betas = fin.read().split('\n')
fin.close()

w = ''
bexec = '/home/ushnish/research/dmcsns/models/SNS/d1b/DEX.branch'
i = 0

nproc = 2
for beta in betas:
	try:
		float(beta)
		sdir = 'set.' + str(i)
		os.system('mkdir ' + sdir)
		os.system('cp -r params ' + sdir + '/')

		w += 'cd ../' + sdir + '\n'
		w += 'mpirun -np 2 ' + bexec + ' ' + './params/params.txt ' + str(beta) + ' >> out &\n\n'
		i += 1
		if((i*nproc)%30 == 0):
			w += 'wait\n'
	except:
		pass		
w += 'wait\n'
print w

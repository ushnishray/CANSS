import os
import sys

fin = open(sys.argv[1])
ints = fin.read().split('\n')
fin.close()

pfile = open(sys.argv[2])
params = pfile.read().split('\n')
pfile.close()

w = ''
#bexec = '/home/ushnish/research/Development/checkouts/DMCsns/cppmods/models/SNS/d1b/DEX.nb'
bexec = '/home/ushnish/research/Development/checkouts/DMCsns/cppmods/models/DBM/d1/DEX.nb'
i = 0

nproc = 2
tproc = 8
for intt in ints:
	try:
		float(intt)
		sdir = 'set.' + str(i)
		os.system('mkdir ' + sdir)
		os.system('mkdir ' + sdir + '/params/')

#		print sdir
		params[int(sys.argv[3])] = str(intt)
		pfile = open(sdir + '/params/param.txt','w')
		for piece in params:
			pfile.write(piece + '\n')
		pfile.close()

		w += 'cd ../' + sdir + '\n'
		w += 'mpirun -np ' + str(nproc) + ' ' + bexec + ' ' + './params/param.txt ' + ' >> out &\n\n'
		i += 1
		if((i*nproc)%tproc == 0):
			w += 'wait\n'
	except:
		pass		
w += 'wait\n'
print w

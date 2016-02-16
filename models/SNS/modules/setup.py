import os

w = ''
bexec = '/home/ushnish/research/Development/checkouts/DMCsns/cppmods/models/SNS/d1b/DEX.branch'
for i in range(0,21):
	beta = -0.1 + 0.01*i
	sdir = 'set.' + str(i)
	os.system('mkdir ' + sdir)
	os.system('cp -r params ' + sdir + '/')
	

	w += 'cd ../' + sdir + '\n'
	w += 'mpirun -np 5 ' + bexec + ' ' + './params/param.txt ' + str(beta) + ' >> out\n\n'
	
print w

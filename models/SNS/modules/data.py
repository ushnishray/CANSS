import os

w = ''
bexec = '/home/ushnish/research/Development/checkouts/DMCsns/cppmods/models/SNS/d1b/DEX.branch'
for i in range(0,21):
	beta = -0.1 + 0.01*i
	sdir = 'set.' + str(i)
#	os.system('mkdir ' + sdir)
#	os.system('cp -r params ' + sdir + '/')
	

#	w += 'cd ../' + sdir + '\n'
#	w += 'mpirun -np 5 ' + bexec + ' ' + './params/param.txt ' + str(beta) + ' >> out\n\n'

	w +=  "gawk '{s+=$3;s2+=$3*$3;n+=1} END {s/=n;s2/=n;se=((s2-s*s)/(n-1))^0.5;print s,se}' " + sdir + '/basicOutput >> Q\n'	
	w +=  "gawk '{s+=$2;s2+=$2*$2;n+=1} END {s/=n;s2/=n;se=((s2-s*s)/(n-1))^0.5;print s,se}' " + sdir + '/basicOutput >> F\n'	
	w +=  "gawk '{s+=$6;s2+=$6*$6;n+=1} END {s/=n;s2/=n;se=((s2-s*s)/(n-1))^0.5;print s,se}' " + sdir + '/basicOutput >> V\n'	

print w

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
		i += 1			
		w +=  "gawk '{s+=$3;s2+=$3*$3;n+=1} END {s/=n;s2/=n;se=((s2-s*s)/(n-1))^0.5;print s,se}' " + sdir + '/basicOutput >> Q\n'	
		w +=  "gawk '{s+=$2;s2+=$2*$2;n+=1} END {s/=n;s2/=n;se=((s2-s*s)/(n-1))^0.5;print s,se}' " + sdir + '/basicOutput >> F\n'	
		w +=  "gawk '{s+=$6;s2+=$6*$6;n+=1} END {s/=n;s2/=n;se=((s2-s*s)/(n-1))^0.5;print s,se}' " + sdir + '/basicOutput >> V\n'	
	except:
		pass		
print w

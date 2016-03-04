import os
import sys

# Use this script to generate Makefile library for every file in the specified directory list
dirfiles = sys.argv[1]
outputfile = sys.argv[2]

dirlist = []
dn = 0
f = open(dirfiles,'r')
for line in f:
	dirpath = line.split()[0]
	dirlist.append(dirpath)
	dn += 1
f.close()

output = 'CPP=mpicxx\n'
output += 'CFLAGS= -O3 -fopenmp -lm -std=gnu++0x -I ./include -I ./core  -I ./movers -I ./waveFunctions -I ./observables -I ./drivers -I ./runners -I ./walkers -I ./Serializer\n'
output += 'LDFLAGS = -fopenmp -lgsl -lgslcblas\n'
output += 'OBJTARGET = ./objLibrary\n'
output += 'OBJECTS=$(OBJTARGET)/*.o\n'
output += 'SLINKER=ar\n'
output += 'STATICLINK=libdmcsns.a\n\n'
output += '#=========================================================================\n\n'

objFileList = ''
rules = ''
for i in range(0,dn):

	dd = dirlist[i]
	
	for subd,d,files in os.walk(dd):
		for ff in files:

			sd = subd[len(dd)+1:]
			
			#Ignore sub-directory that contains test

			if(ff[-4:] == '.cpp' and not('test' in sd)):
				print dd + '/' + ff[:-4] + '.o',sd
				objfile = ff[:-4] + '.o'
								
				cppff = dd + '/' + ff		
				hff = dd + '/' + ff[:-4] + '.h'
				if(~os.path.isfile(hff)):
					hff = ''
				off = '$(OBJTARGET)/' + objfile

				objFileList += off + ' '

				rules += off + ' : ' + cppff + ' ' + hff + '\n'
				rules += '\t$(CPP) $(CFLAGS) -c $(DEBUG) ' + cppff + '\n'
				rules += '\tmv ' + objfile + ' $(OBJTARGET)/\n\n' 

output += '$(STATICLINK): ' + objFileList + '\n'
output += '\t$(SLINKER) -rcs $(STATICLINK) $(OBJECTS)\n\n'

output += 'clean: \n'
output += '\t rm -f $(OBJECTS)\n'
output += '\t rm $(STATICLINK)\n'
output += '#=========================================================================\n\n'
output += rules

f = open(outputfile,'w')
f.write(output)
f.close()

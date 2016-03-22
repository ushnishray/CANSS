BASEPATH=../../../
SSEPATH=-I $(BASEPATH)/include -I $(BASEPATH)/core -I $(BASEPATH)/movers -I $(BASEPATH)/observables -I $(BASEPATH)/drivers -I $(BASEPATH)/waveFunctions -I$(BASEPATH)/walkers -I$(BASEPATH)/runners -I$(BASEPATH)/Serializer
CPP=mpicxx 
FORTRAN=
FLAGS=
CFLAGS= -O2 -fopenmp -lm -std=gnu++0x -lgsl $(SSEPATH) -I ./ 
LDFLAGS = -fopenmp -lgsl -lgslcblas -L$(BASEPATH) -ldmcsns 
F_SOURCES= 
DEBUG= 
CPP_SOURCES= ./*.cpp  
OBJECTS=*.o
EXECUTABLE =  DEX.branch

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJECTS): $(CPP_SOURCES)
	$(CPP) $(CFLAGS) -c $(DEBUG) $(CPP_SOURCES)
#	$(FORTRAN) $(FFLAGS) -c $(F_SOURCES)

clean: 
	rm -f *.o

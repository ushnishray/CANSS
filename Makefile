CPP=mpicxx
CFLAGS= -O3 -fopenmp -lm -std=gnu++0x -I ./include -I ./core  -I ./movers -I ./waveFunctions -I ./observables -I ./drivers -I ./runners -I ./walkers -I ./Serializer
LDFLAGS = -fopenmp -lgsl -lgslcblas
OBJTARGET = ./objLibrary
OBJECTS=$(OBJTARGET)/*.o
SLINKER=ar
STATICLINK=libdmcsns.a

#=========================================================================

$(STATICLINK): $(OBJTARGET)/Walker.o $(OBJTARGET)/Qhistogram.o $(OBJTARGET)/AutoCorrI.o $(OBJTARGET)/Whistogram.o $(OBJTARGET)/Phistogram.o $(OBJTARGET)/AutoCorrFull.o $(OBJTARGET)/BasicObs.o $(OBJTARGET)/CloneMultiplicity.o $(OBJTARGET)/Density.o $(OBJTARGET)/AutoCorr.o $(OBJTARGET)/branchMechanismLimited.o $(OBJTARGET)/branchLocal.o $(OBJTARGET)/runstratwb.o $(OBJTARGET)/runstratnb.o $(OBJTARGET)/MPIBasicRunner.o $(OBJTARGET)/branchMechanismFull.o 
	$(SLINKER) -rcs $(STATICLINK) $(OBJECTS)

clean: 
	 rm -f $(OBJECTS)
	 rm $(STATICLINK)
#=========================================================================

$(OBJTARGET)/Walker.o : core/Walker.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) core/Walker.cpp
	mv Walker.o $(OBJTARGET)/

$(OBJTARGET)/Qhistogram.o : observables/Qhistogram.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/Qhistogram.cpp
	mv Qhistogram.o $(OBJTARGET)/

$(OBJTARGET)/AutoCorrI.o : observables/AutoCorrI.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/AutoCorrI.cpp
	mv AutoCorrI.o $(OBJTARGET)/

$(OBJTARGET)/Whistogram.o : observables/Whistogram.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/Whistogram.cpp
	mv Whistogram.o $(OBJTARGET)/

$(OBJTARGET)/Phistogram.o : observables/Phistogram.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/Phistogram.cpp
	mv Phistogram.o $(OBJTARGET)/

$(OBJTARGET)/AutoCorrFull.o : observables/AutoCorrFull.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/AutoCorrFull.cpp
	mv AutoCorrFull.o $(OBJTARGET)/

$(OBJTARGET)/BasicObs.o : observables/BasicObs.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/BasicObs.cpp
	mv BasicObs.o $(OBJTARGET)/

$(OBJTARGET)/CloneMultiplicity.o : observables/CloneMultiplicity.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/CloneMultiplicity.cpp
	mv CloneMultiplicity.o $(OBJTARGET)/

$(OBJTARGET)/Density.o : observables/Density.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/Density.cpp
	mv Density.o $(OBJTARGET)/

$(OBJTARGET)/AutoCorr.o : observables/AutoCorr.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) observables/AutoCorr.cpp
	mv AutoCorr.o $(OBJTARGET)/

$(OBJTARGET)/branchMechanismLimited.o : runners/branchMechanismLimited.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) runners/branchMechanismLimited.cpp
	mv branchMechanismLimited.o $(OBJTARGET)/

$(OBJTARGET)/branchLocal.o : runners/branchLocal.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) runners/branchLocal.cpp
	mv branchLocal.o $(OBJTARGET)/

$(OBJTARGET)/runstratwb.o : runners/runstratwb.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) runners/runstratwb.cpp
	mv runstratwb.o $(OBJTARGET)/

$(OBJTARGET)/runstratnb.o : runners/runstratnb.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) runners/runstratnb.cpp
	mv runstratnb.o $(OBJTARGET)/

$(OBJTARGET)/MPIBasicRunner.o : runners/MPIBasicRunner.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) runners/MPIBasicRunner.cpp
	mv MPIBasicRunner.o $(OBJTARGET)/

$(OBJTARGET)/branchMechanismFull.o : runners/branchMechanismFull.cpp 
	$(CPP) $(CFLAGS) -c $(DEBUG) runners/branchMechanismFull.cpp
	mv branchMechanismFull.o $(OBJTARGET)/



INITIAL  = mmsn
HYDRO    = euler
BOUNDARY = fixed
OUTPUT   = h5out
RESTART  = h5in
PLANETS  = circular

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /home/install/app/hdf5
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/local
endif

CC = mpicc
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lm -lhdf5

OBJ = main.o readpar.o timestep.o onestep.o riemann.o mpisetup.o gridsetup.o domain.o misc.o geometry.o faces.o exchange.o plm.o report.o profiler.o planet.o omega.o analysis.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o $(BOUNDARY).o $(RESTART).o $(PLANETS).o #snapshot.o

default: disco

%.o: %.c paul.h
	$(CC) $(FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(PLANETS).o : Planets/$(PLANETS).c paul.h
	$(CC) $(FLAGS) $(INC) -c Planets/$(PLANETS).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c paul.h
	$(CC) $(FLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c paul.h
	$(CC) $(FLAGS) $(INC) -c Restart/$(RESTART).c

disco: $(OBJ) paul.h
	$(CC) $(FLAGS) $(LIB) -o disco $(OBJ)

clean:
	rm -f *.o disco

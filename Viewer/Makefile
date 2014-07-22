
EXE  = vdisco

GLUT = /usr/X11R6
HDF5 = /opt/local

LIBS = -L$(GLUT)/lib -L$(HDF5)/lib -framework GLUT -framework OpenGL -lXmu -lXext -lXi -lX11 -lm -lhdf5 -lz
INCL = -I$(GLUT)/include -I$(HDF5)/include

CC = g++

default: $(EXE)

$(EXE): $(EXE).c
	$(CC) -c colorbar.c
	$(CC) -O3 -Wall -o $@ $< colorbar.o $(INCL) $(LIBS)

clean:
	rm -f *.o $(EXE)

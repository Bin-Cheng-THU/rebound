include Makefile.defs

OPT+= -fPIC -DLIBREBOUND

SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_leapfrog.c boundary.c input.c output.c collision.c communication_mpi.c display.c tools.c simulationarchive.c glad.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.so 

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) -o $@ $<

librebound.so: $(OBJECTS)
	@echo ""        
	@echo "Linking share library $@ ..."
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@echo "Removing object files *.o ..."
	@-rm -f *.o
	@echo "Removing shared library librebound.so ..."
	@-rm -f *.so
	@echo "Removing coverage files ..."
	@-rm -f *.gcda
	@-rm -f *.gcno
	

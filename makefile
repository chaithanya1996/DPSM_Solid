#CC = g++
CC = clang++
CFLAGS= -larmadillo  -g -fopenmp -O3
#WORKING_DIRECTORY = $(pwd)

all: DPSM.cpp helper.o core.o interface.o 
	#$(CC) DPSM.cpp dpsm_helpers.o  dpsm_core.o dpsm_interface.o -o DPSM.o  $(CFLAGS)
helper.o: dpsm_helpers.cpp dpsm_helpers.hpp
	$(CC) -c dpsm_helpers.cpp -o dpsm_helpers.o $(CFLAGS)
core.o: dpsm_core.cpp dpsm_core.hpp
	$(CC) -c dpsm_core.cpp dpsm_helpers.o -o dpsm_core.o $(CFLAGS)
interface.o: dpsm_interface.cpp dpsm_interface.hpp
	$(CC) -c dpsm_interface.cpp dpsm_helpers.o dpsm_core.o -o dpsm_interface.o $(CFLAGS)
clean:
	rm *.o


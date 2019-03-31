CC = g++
#CC = clang++
INTEL_MKL= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl 
CFLAGS=  -lopenblas -llapack  -larmadillo  -g -fopenmp -O3 #$(INTEL_MKL) -m64 -I${MKLROOT}/include
#WORKING_DIRECTORY = $(pwd)

all: DPSM.cpp helper.o core.o interface.o 
	$(CC) DPSM.cpp dpsm_helpers.o  dpsm_core.o dpsm_interface.o -o DPSM.o  $(CFLAGS)
helper.o: dpsm_helpers.cpp dpsm_helpers.hpp
	$(CC) -c dpsm_helpers.cpp -o dpsm_helpers.o $(CFLAGS)
core.o: dpsm_core.cpp dpsm_core.hpp
	$(CC) -c dpsm_core.cpp dpsm_helpers.o -o dpsm_core.o $(CFLAGS)
interface.o: dpsm_interface.cpp dpsm_interface.hpp
	$(CC) -c dpsm_interface.cpp dpsm_helpers.o dpsm_core.o -o dpsm_interface.o $(CFLAGS)
clean:
	rm *.o
	rm *.csv

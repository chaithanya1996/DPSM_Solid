CC = g++
#CC = clang++
CFLAGS=  -lopenblas -llapack  -larmadillo  -g -fopenmp -O3 #$(INTEL_MKL) -m64 -I${MKLROOT}/include
#WORKING_DIRECTORY = $(pwd)

all: DPSM.cpp helper.o core.o 
	$(CC) DPSM.cpp dpsm_helpers.o  dpsm_core.o  -o DPSM.o  $(CFLAGS)
EXP: DPSM_time.cpp helper.o core.o 
	$(CC) DPSM_time.cpp dpsm_helpers.o  dpsm_core.o  -o DPSM_time.o  $(CFLAGS)
TEST: DPSM_testing.cpp helper.o core.o
	$(CC) DPSM_testing.cpp dpsm_helpers.o  dpsm_core.o    -o DPSM_testing.o  $(CFLAGS)
POST: DPSM_time_post.cpp helper.o core.o 
	$(CC) DPSM_time_post.cpp dpsm_helpers.o  dpsm_core.o    -o DPSM_time_post.o  $(CFLAGS)
helper.o: dpsm_helpers.cpp dpsm_helpers.hpp
	$(CC) -c dpsm_helpers.cpp -o dpsm_helpers.o $(CFLAGS)
core.o: dpsm_core.cpp dpsm_core.hpp
	$(CC) -c dpsm_core.cpp dpsm_helpers.o -o dpsm_core.o $(CFLAGS)
# interface.o: dpsm_interface.cpp dpsm_interface.hpp
# 	$(CC) -c dpsm_interface.cpp dpsm_helpers.o dpsm_core.o -o dpsm_interface.o $(CFLAGS)
clean:
	rm *.o
	rm *.csv

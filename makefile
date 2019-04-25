CC = g++
#CC = clang++
CFLAGS=  -lopenblas -llapack  -larmadillo  -g -fopenmp -O3 #$(INTEL_MKL) -m64 -I${MKLROOT}/include
#WORKING_DIRECTORY = $(pwd)

all: DPSM.cpp
	$(CC) DPSM.cpp   -o DPSM.o  $(CFLAGS)
EXP: DPSM_time.cpp 
	$(CC) DPSM_time.cpp  -o DPSM_time.o  $(CFLAGS)
TEST: DPSM_testing.cpp 
	$(CC) DPSM_testing.cpp   -o DPSM_testing.o  $(CFLAGS)
#POST: DPSM_time_post.cpp
#	$(CC) DPSM_time_post.cpp    -o DPSM_time_post.o  $(CFLAGS)
POST:DPSM_post.cpp
	$(CC) DPSM_post.cpp -o DPSM_post.o $(CFLAGS)
clean:
	rm *.o
	rm *.csv

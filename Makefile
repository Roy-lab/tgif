#CHANGE PATHS AS NEEDED:
INCLUDE_PATH = ./gsl/include
LIBRARY_PATH = ./gsl/lib
#INCLUDE_PATH = ${CONDA_PREFIX}/include
#LIBRARY_PATH = ${CONDA_PREFIX}/lib

#compiler: gcc for C programs, g++ for C++ programs
XX = g++
CC = gcc

#compiler flags
CFLAGS = -g
GSLFLAGS = -lgsl -lgslcblas

all: clean tgif 

matf:
	$(CC) -c -o modules/random_svd/matrix_funcs.o modules/random_svd/matrix_vector_functions_gsl.c -I${INCLUDE_PATH}

rsvd:
	$(CC) -c -o modules/random_svd/rsvd.o modules/random_svd/low_rank_svd_algorithms_gsl.c -I${INCLUDE_PATH}

tgif:
	$(XX) tgif-db.cpp modules/*.cpp modules/random_svd/*.o -o tgif-db $(CFLAGS) -L${LIBRARY_PATH} ${GSLFLAGS} -I${INCLUDE_PATH}
	$(XX) tgif-dc.cpp modules/*.cpp modules/random_svd/*.o -o tgif-dc $(CFLAGS) -L${LIBRARY_PATH} ${GSLFLAGS} -I${INCLUDE_PATH}

clean:
	rm tgif-db
	rm tgif-dc

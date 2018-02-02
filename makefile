#LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:/usr/lib64/:/usr/lib64/psm2
#export LD_LIBRARY_PATH

assn2: powermethod.c functions.c powermethod.h
	mpicc -o assn2 powermethod.c functions.c -lm

clean: 
	-rm *.o
	-rm assn2

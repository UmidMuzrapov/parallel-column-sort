## Shared-Memory Parallel Column Sort
The program implements the parallelized version of a column sorts that sorts nonnegative integers using pthreads. 

## Requirements
C99 and pthread.h

## Installation
Include columSort.h in the user program.
```
#include "columnSort.h"
```
Compile and run your program. This is an example Makefile:
```
parsort:	threadColumnSort.o columnSortHelper.o  main.o
	gcc -o parsort threadColumnSort.o columnSortHelper.o main.o -lm -lpthread
	
main.o:	main.c columnSort.h columnSortHelper.h
	gcc -c -O2 -std=c99 main.c

columnSortHelper.o:	columnSortHelper.c columnSortHelper.h 
	gcc -c -O2 -std=c99 columnSortHelper.c

threadColumnSort.o:	threadColumnSort.c columnSortHelper.h
	gcc -c -O2 -std=c99 threadColumnSort.c

clean:
	rm -f *.o parsort
```

## API
```
void columnSort(int *A, int numThreads, int length, int width, double *elapsedTime)
```
Input
A: input array to sort
numThreads: number of light-weight processes to use
length: number of rows
width: number of columns
elapsedTime: used to pass back the elapsed time to do column-sort.

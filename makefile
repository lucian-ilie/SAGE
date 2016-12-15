SHELL=/bin/sh

#Compilers
CC=gcc 
CC2=g++

#FLAGS
MYFLAGS = -O3 -DNDEBUG
CFLAGS = -W -Wall -Winline -O3

.PHONY: all
#Pattern rule for all objects files
#%.o: %.c *.h
#	$(CC) -c $(CFLAGS) $< -o $@

#For making all object file and executable files
SAGE: sage.o cs2.o smithWaterman.o kseq.h
	$(CC2) $(CFLAGS) -o SAGE sage.o cs2.o smithWaterman.o -lm -lz

#For cs2
cs2.o: CS2/cs2.c CS2/parser_cs2.c CS2/types_cs2.h CS2/timer.c 
	$(CC) $(MYFLAGS) -c CS2/cs2.c

smithWaterman.o : smithWaterman.c
	$(CC) $(CFLAGS) -c smithWaterman.c
  
sage.o : sage.c
	$(CC) $(CFLAGS) -c sage.c

#For removing all object files and executable files
clean: 
	rm -f *.o *.a SAGE





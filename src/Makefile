CC=	gcc
CFLAG= -O3
SRCS=	util_lib.c hmm_lib.c run_hmm.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o

SRCSPOOL=   util_lib.c hmm_lib.c run_hmm_pool.c thpool.c
OBJPOOL=    util_lib.o hmm_lib.o run_hmm_pool.o thpool.o

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC) $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o -lm -lpthread -Wall -pg -pedantic-errors

pool:   $(OBJPOOL)
	$(CC) $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm_pool.o thpool.o -lm -lpthread

sanitize:	$(OBJ)
	$(CC) $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o -lm -lpthread -Wall -pg -fsanitize=address

clean:
	rm -rf example/*.tmp.*
	rm -rf *.o FragGeneScan* *~

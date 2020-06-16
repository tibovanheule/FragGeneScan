CC=	gcc
CFLAG= -Ofast
SRCS=	util_lib.c hmm_lib.c run_hmm.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o

SRCSPOOL=   util_lib.c hmm_lib.c run_hmm_pool.c thpool.c
OBJPOOL=    util_lib.o hmm_lib.o run_hmm_pool.o thpool.o

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC) $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o -lm -lpthread

pool:   $(OBJPOOL)
	$(CC) $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm_pool.o thpool.o -lm -lpthread

clean:
	rm -rf example/*.tmp.*
	rm -rf *.o FragGeneScan* *~

CC=	gcc
CFLAG= -Ofast
SRCS=	util_lib.c hmm_lib.c run_hmm.c sprint_irc.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o sprintf_irc.o

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC) $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o sprintf_irc.o  -lm -lpthread -fsanitize=address

clean:
	rm -rf *.o FragGeneScan* *~

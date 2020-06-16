CC=	gcc
CFLAG= -O3
SRCS=	util_lib.c
OBJ=	util_lib.o 

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC)  $(CFLAG) -o MainSwitchTest util_lib.o -lm -lpthread

clean:
	rm -rf *.o FragGeneScan* *~

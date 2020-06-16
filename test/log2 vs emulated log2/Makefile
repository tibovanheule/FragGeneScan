CC=	gcc
CFLAG= -O3
SRCS=	main.c
OBJ=	main.o 

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC)  $(CFLAG) -o Mainlog2Test main.o -lm -lpthread

clean:
	rm -rf *.o Mainlog2Test* *~

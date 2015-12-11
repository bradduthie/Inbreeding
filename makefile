CXX    = icc
CFLAGS = -fopenmp -ansi -Wall -pedantic -std=c99

RunInbr:	InMult.o Inbreed.o retain.o mortality.o parents.o babies.o mateselect.o 		deadbottom.o Rbuild.o vectorrank.o arraysort2D.o landscape.o randnormINT.o 			sumstats.o randnorm.o randunif.o FFT.o as183.o inbrdep.o
	$(CC) $(CFLAGS) -o InMult InMult.o Inbreed.o retain.o mortality.o parents.o 			babies.o mateselect.o deadbottom.o Rbuild.o vectorrank.o arraysort2D.o 			landscape.o randnormINT.o sumstats.o randnorm.o randunif.o FFT.o as183.o 			inbrdep.o -lm

InMult.o: InMult.c
	$(CC) $(CFLAGS) -c InMult.c

Inbreed.o: Inbreed.c
	$(CC) $(CFLAGS) -c Inbreed.c

retain.o: retain.c
	$(CC) $(CFLAGS) -c retain.c

mortality.o: mortality.c
	$(CC) $(CFLAGS) -c mortality.c

parents.o: parents.c
	$(CC) $(CFLAGS) -c parents.c

babies.o: babies.c
	$(CC) $(CFLAGS) -c babies.c

mateselect.o: mateselect.c
	$(CC) $(CFLAGS) -c mateselect.c

deadbottom.o: deadbottom.c
	$(CC) $(CFLAGS) -c deadbottom.c

Rbuild.o: Rbuild.c
	$(CC) $(CFLAGS) -c Rbuild.c

vectorrank.o: vectorrank.c
	$(CC) $(CFLAGS) -c vectorrank.c

arraysort2D.o: arraysort2D.c
	$(CC) $(CFLAGS) -c arraysort2D.c

landscape.o: landscape.c
	$(CC) $(CFLAGS) -c landscape.c

randnormINT.o: randnormINT.c
	$(CC) $(CFLAGS) -c randnormINT.c

sumstats.o: sumstats.c
	$(CC) $(CFLAGS) -c sumstats.c

randnorm.o: randnorm.c
	$(CC) $(CFLAGS) -c randnorm.c

randunif.o: randunif.c
	$(CC) $(CFLAGS) -c randunif.c

FFT.o: FFT.c
	$(CC) $(CFLAGS) -c FFT.c

as183.o: as183.c
	$(CC) $(CFLAGS) -c as183.c

inbrdep.o: inbrdep.c
	$(CC) $(CFLAGS) -c inbrdep.c

clean:
	rm -rf *.o *.txt daj
	find . -type f -name out\* -exec rm {} \;

VERSION=10.4
MATHLINKPATH = /usr/local/Wolfram/Mathematica/${VERSION}/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
SYS = Linux

MPREP = $(MATHLINKPATH)/mprep

CC = /usr/bin/gcc

LIBFLAGS=-L$(MATHLINKPATH)
FLAGS=-I$(MATHLINKPATH) -O2 -Wall -std=c99 -fno-common 
EXTRA_CFLAGS=
OMP_FLAGS=-fopenmp
#OMP_FLAGS=

OBJECTS = GaussTerminal.o GaussInitData.o GaussCommon.o GaussUserProblem.o GaussCoefficients.o 
REALOBJECTS = math-Gausstm.o math-Gauss.o GaussCommon.o GaussUserProblem.o GaussCoefficients.o 
QUADOBJECTS = quad-Gausstm.o quad-Gauss.o GaussCommon.o GaussUserProblem.o GaussCoefficients.o 
OBJECTS3 = math-GaussParatm.o math-GaussPar.o GaussCommon.o GaussUserProblem.o  GaussCoefficients.o 

all : 
	make clean;
	make term-Gauss;
	make math-Gauss;
	make quad-Gauss;
	make math-GaussPar

term-Gauss: $(OBJECTS) 
	${CC} $(OBJECTS) ${OMP_FLAGS} -lm -lquadmath -o GaussTerminal.exe

math-Gausstm.c : math-Gauss.tm
	${MPREP} math-Gauss.tm -o math-Gausstm.c

quad-Gausstm.c : quad-Gauss.tm
	${MPREP} quad-Gauss.tm -o quad-Gausstm.c

math-GaussParatm.c : math-GaussPar.tm
	${MPREP} math-GaussPar.tm -o math-GaussParatm.c

%.o: %.c prec.h def.h 
	${CC} ${FLAGS} ${OMP_FLAGS} -c -I./  $<

math-Gauss : $(REALOBJECTS)    
	${CC}  $(REALOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -o $@ 

quad-Gauss : $(QUADOBJECTS)   
	${CC} $(QUADOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -lquadmath -o $@ 

math-GaussPar : $(OBJECTS3) 
	${CC} $(OBJECTS3) $(LIBFLAGS) ${OMP_FLAGS}  -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -lquadmath -o $@

clean :
	rm -rf *.o *~ core *tm.* *.exe 





